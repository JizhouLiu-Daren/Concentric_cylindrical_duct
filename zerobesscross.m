function k = zerobesscross(nu,alpha,N,varargin)
    % k = zerocross(m,alpha,n)
    % Returns the first N roots of the Bessel function cross products of order m.
    %
    % Optional arguments:
    %  - Boundary conditions: one of
    %    'DD' (default):
    %      J(m,x)*Y(m,alpha*x) - Y(m,x)*J(m,alpha*x) = 0,
    %    'DN':
    %      J(m,x)*Y'(m,alpha*x) - Y(m,x)*J'(m,alpha*x) = 0,
    %    'ND':
    %      J(m,alpha*x)*Y'(m,x) - Y(m,alpha*x)*J'(m,x) = 0,
    %    'NN':
    %      J'(m,x)*Y'(m,alpha*x) - Y'(m,x)*J'(m,alpha*x) = 0,
    %
    %  - Relative tolerance (default 10*eps)
    
    numvarargs = length(varargin);
    optargs = {'DD', 10*eps}; %Defaults
    optargs(1:numvarargs) = varargin;
    [T,tol] = optargs{:};

    if alpha < 1
        a = 1/alpha;
        if strcmp(T, 'DN')
            T = 'ND';
        elseif strcmp(T,'ND')
            T = 'DN';
        end
    else
        a = alpha;
    end


    J = @besselj;
    Y = @bessely;
    dJ = @(mu,x) (J(mu-1,x) - J(mu+1,x))/2;
    dY = @(mu,x) (Y(mu-1,x) - Y(mu+1,x))/2;

    %True phase function
    function [t,dt] = th(x)
        j = J(nu,x);
        y = Y(nu,x);

        %Fix matlab bug.
        y(~isfinite(y)) = -Inf;
        
        t = atan2(y,j);
        dt = 2./(pi*x*(j.^2+y.^2));
    end

    jp2 = nu + 2.57809*nu^(1/3); % lower bound for 2nd root of J'
    %Approximate phase function and its derivative.
    %From Olver's expansion
    function [t,dt] = tha(x)
        if x >= jp2
            r = sqrt(x^2-nu^2);
            xi = (r-nu*asec(x/nu));
            t = xi-pi/4; 
            dt = r/x;
        else
            [t,dt] = th(x);
        end
    end

    %True phase function for derivatives.
    function [t,dt] = ph(x)
        dj = dJ(nu,x);
        dy = dY(nu,x);

        %Fix a matlab bug as well as our evaluation method 
        dy(~isfinite(dy)) = Inf; 
        
        t = atan2(dy,dj);
        dt = 2*(x.^2 - nu^2)./(pi*x.^3.*(dj.^2+dy.^2)); 
    end

    
    %approximation for y'_1 taken from zerobess.
    yp1 = nu + 1.8212*(nu+1)^(1/3)  + 0.9224*(nu+1)^(-1/3) - 0.5561*(nu+1)^(-2/3);

    %Approximate derivative phase function.
    %From Olver's expansion
    function [t,dt] = pha(x)
        if x >= yp1
            r = sqrt(x^2-nu^2);
            xi = (r-nu*asec(x/nu));
            dt = r/x;
            t = xi + pi/4;
        else
            [t,dt] = ph(x);
        end
    end

    k = zeros(N,1);

    l0=0;

    if T(1) == 'D'
        ph2 = @th;
        pha2 = @tha;
    elseif T(1) == 'N'
        ph2 = @ph;
        pha2 = @pha;
        % If the first phase type is derivative, first zeros occurs when total phase = 0;
        l0 = -1; 
    else
        error('Type must be either ''D'' or ''N''.') 
    end

    if T(2) == 'D'
        ph1 = @th;
        pha1 = @tha;
    elseif T(2) == 'N'
        ph1 = @ph;
        pha1 = @pha;
    else
        error('Type must be either ''D'' or ''N''.') 
    end

    %Extra root vanishes at nu==0 for NN case
    if T(1) == 'N' && T(2) == 'N' && nu == 0
        l0 = 0;
    end


    for n = 1:N
        %Initial guess 
        m = n+l0;
        k0 = m*pi/(a-1);

        % Guess above is a bad start for small roots
        if k0 <= nu/a
           %Very rough guess
           k0 = (nu + (m+1)*(nu+1)^(1/3))/a;
        end


        %In the special case 
        wdt = eps(jp2);
        if strcmp(T,'ND') && (tha(a*(yp1-wdt))-pha(yp1-wdt) < m*pi) && (tha(a*(yp1+wdt))-pha(yp1+wdt) > m*pi)
            k0 = yp1;
        elseif strcmp(T,'ND') && (tha(jp2-wdt)-pha(jp2/a-wdt) < m*pi && tha(jp2+wdt)-pha(jp2/a+wdt) > m*pi)
            warning('special case 2')
            k0 = jp2/a;
        else
            err = 1;
            i=0;
            while err > min(1e-4,0.1*pi/a) && i < 4;
                [t1,dt1] = pha1(a*k0);
                [t2,dt2] = pha2(k0);
                update = (t1-t2-m*pi)/(dt1*a-dt2);
                k0 = k0 - update; 
                err = abs(update/k0);
                i=i+1;
            end
            
        end

        err = 1;
        i=0;
        while err > tol && i < 20
            [t1,dt1] = ph1(a*k0);
            [t2,dt2] = ph2(k0);
            update = tan(t1-t2)./(a*dt1-dt2);
            k0 = k0 - update;
            err = abs(update/k0);
            i=i+1;
        end
        if i == 20
            warning('unable to reach tol at n = %d',n)
        end
        k(n) = k0;
    end

    %Transform roots if nessecary
    if alpha < 1
        k = k/alpha;
    end

end

