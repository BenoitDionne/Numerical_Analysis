% X1 = cobweb(funct, funct_name, a, b, I, VX0, N, UNIT, V, index, W)
%
% Program to draw the iteration of a discrete dynamical system.
% It also returns the last value of the iteration.
%
% Given:
%  funct: the function handle that define the discrete dynamical system
%     x_{n+1} = f(x_n , a, b)
%  funct_name: the title for the graph.
%  a , b:  the parameters allowed in the function.
%  I:  the interval such that  f:I -> I .
%  VX0: the initial conditions.
%  N: the number of iterations.
%  UNIT: the number of points per unit length.  Increase this number
%     to get nicer curves.  The default is 100.
%  V: the name of the dependant variable.  The default is x .
%  index: the name of the index.  The default is n .
%  W: width of the arrow.  Set to 0 to not draw the arrow.
%
function X1 = cobweb(funct, funct_name, a, b, I, VX0, N, UNIT, V, index, W)

    % default arguments
    arguments
        funct;
        funct_name string;
        a double;
        b double;
        I (1,2) double;
        VX0 (1,:) double;
        N int32;
        UNIT int32 = 100;
        V string = 'x';
        index string = 'n';
        W double = -1;
    end

    % Set up the graph
    cla;
    axis equal;
    grid on;
    hold on;
    xlabel( strcat(V,'_',index) );
    ylabel( strcat(V,'_{',index,'+1}') );
    title(funct_name);
    figure(gcf);

    % Number of initial conditions
    M = size(VX0,2);

    % To draw the arrows
    S0 = zeros(1,M);
    S1 = zeros(1,M);

    % Find the largest and smallest coordinates.
    MAX = VX0(1,1);
    MIN = VX0(1,1);

    for m=1:M
        % initial condition
        X0 = VX0(1,m);
        if ( X0 > MAX ) MAX = X0; end
        if ( X0 < MIN ) MIN = X0; end

        % Save X0 for the end to draw the arrow
        S0(1,m) = X0;

        % First iteration
        X1 = feval(funct,X0,a,b);
        if ( X1 > MAX ) MAX = X1; end
        if ( X1 < MIN ) MIN = X1; end

        % Save X1 for the end to draw the arrow
        S1(1,m) = X1;

        y = linspace(0, X1, 2);
        x = ones(2,1)*X0;
        plot(x,y,'r');
        x = linspace(X0, X1, 2);
        y = ones(2,1)*X1;
        plot(x,y,'r');
        X0 = X1;
        for i=1:N
            % Follwing iterations
            X1 = feval(funct,X0,a,b);
            if ( X1 > MAX ) MAX = X1; end
            if ( X1 < MIN ) MIN = X1; end
            y = linspace(X0, X1, 2);
            x = ones(2,1)*X0;
            plot(x,y,'r');
            x = linspace(X0, X1, 2);
            y = ones(2,1)*X1;
            plot(x,y,'r');
            X0 = X1;
        end
    end

    % plot the arrow(s)
    if ( W < 0 )
        W = 0.02*(MAX - MIN);
    end;
  
    for m=1:M
        s = [S0(1,m) - W, S0(1,m)];
        if ( S1(1,m) < 0 )
            t = [S1(1,m)/2 + W, S1(1,m)/2];
        else
            t = [S1(1,m)/2 - W, S1(1,m)/2];
        end
        plot(s,t,'r')
        s = [S0(1,m) + W, S0(1,m)];
        plot(s,t,'r')
    end

    % plot the x and y axes
    if ( I(2) > MAX ) MAX = I(2); end
    if ( I(1) < MIN ) MIN = I(1); end
    if ( MIN <= 0 && MAX >= 0 )
        x = linspace(MIN, MAX, 2);
        y = 0*x;
        plot(x,y,'k');
        plot(y,x,'k');
    end

    % plot the line y = x
    x = linspace( MIN, MAX, 2);
    plot(x,x,'b');

    % plot the graph of the function y = f(x,a,b)
    x = linspace( MIN, MAX, numb_points( MAX-MIN , UNIT) );
    y = feval(funct,x,a,b);
    plot(x,y,'k');

    % Visible part of the graph
    axis( [ MIN, MAX, MIN, MAX ]);

end  % for the function

% Useful sub-function
function y = numb_points(x , UNIT)
    y = ceil( abs(x)*UNIT );
    if ( y < 10 )
        y = 10;
    end
end   %  for the function
