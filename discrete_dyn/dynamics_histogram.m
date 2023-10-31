%  dynamics_histogram(funct, funct_name, a, b, I, VX0, N, UNIT, V, index, W)
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
%  X0: the initial condition.
%  N: the number of iterations.
%  Nparts: the number of partitions of the interval I.  The default is 20.
%  V: the name of the dependant variable.  The default is x .
%
% Example:
% ff = @(x,a,b) a*x.*(1-x)
% dynamics_histogram(ff,'logistic', 3.5, 0, [0,1],0.5, 100, 20, 'x')
%
function dynamics_histogram(funct, funct_name, a, b, I, X0, N, Nparts, V)

    % default arguments
    arguments
        funct;
        funct_name string;
        a double;
        b double;
        I (1,2) double;
        X0 double;
        N int32;
        Nparts int32 = 20;
        V string = 'x';
     end

    % Set up the histogram
    cla;

    Xi = [funct(X0,a,b)];
    for i=1:1:N-1
        X = funct(Xi(i),a,b);
        if ( X >= I(1) && X <= I(2) )
            Xi = [Xi,X];
        end
    end
    histogram(Xi,Nparts,'Normalization','probability','FaceColor','b','EdgeColor','b');
    xlabel(V);
    ylabel('Percentage');
    title(funct_name);
    grid on;
    axis( [ I(1), I(2), 0, 1 ]);
end

% Using the function bar()
%
% dynamics_histogram(funct, funct_name, a, b, I, VX0, N, UNIT, V, index, W)
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
%  X0: the initial condition.
%  N: the number of iterations.
%  Nparts: the number of partitions of the interval I.  The default is 20.
%  V: the name of the dependant variable.  The default is x .
%
% Example:
% ff = @(x,a,b) a*x.*(1-x)
% dynamics_histogram(ff,'logistic', 3.5, 0, [0,1],0.5, 100, 20, 'x')
%
% function dynamics_histogram(funct, funct_name, a, b, I, X0, N, Nparts, V)
% 
%     % default arguments
%     arguments
%         funct;
%         funct_name string;
%         a double;
%         b double;
%         I (1,2) double;
%         X0 double;
%         N int32;
%         Nparts int32 = 20;
%         V string = 'x';
%     end
%
%     % Set up the histogram
%     cla;
% 
%     data = zeros(1,Nparts);
%     Npartsp1 = Nparts+1;
%     ends = linspace(I(1),I(2),Npartsp1);
%     middles = ends(2:1:Npartsp1) - (ends(2)-ends(1))/2;
%     count = 0;
%
%     for i=1:1:N
%         X0 = funct(X0,a,b);
%         if ( X0 >= I(1) && X0 <= I(2) )
%             count = count + 1;
%             for j=2:1:Npartsp1
%                 if ( ends(j) > X0 )
%                     data(j-1) = data(j-1) + 1
%                     break;
%                 end
%             end
%         end
%     end
%
%     % plot the histogram
%     data = data/count;
%     bar(middles, data, 1);
%     xlabel(V);
%     ylabel('Percentage');
%     title(funct_name);
%    grid on;
%    axis( [ I(1), I(2), 0, 1 ]);
% end
