%  ospline(P,spacing,handle)
%
%  This program plots a spline that shadows the points
%    P(1,:), P(2,:), ... , P(npoints,:).  The o in ospline if
%    for other splines.
%
%  We set PLL = PL = P(1,:).  We begin the plotting with the B-spline
%  associated to  PLL, PL, P(1,:) and P(2,:).  The point P(1,:) will
%  then be on the B-spline for P(i,:) with i=1 to npoints
%  We set PRR = PR = P(npoints,:).  We end the plotting with the
%  B-spline associated to P(npoints-1,:), P(npoints,:), PR and PRR.
%  The point P(npoints,:) will then be on the B-spline for P(i,:) with
%  i=1 to npoints
%
%  The following information must be given by the user:
%    The points P(i,:) for i=1 to npoints,
%    The handle of a figure.
%      If no handle of a figure is given, the program look for
%      an already open figure.  If there is an open figure, the
%      program uses its handle.  If there is no open figure, the
%      program open a new figure.
%    The number N of subintervals of [0,1] for each piecewise
%      curves.  The default is 20.
%
function ospline(P,N,handle)
    % default arguments
    arguments
        P (:,2) double;
        N = 20;
        handle = gcf;
    end

    figure(handle);
    hold on;

    npoints = size(P,1);
    Q = [ P(1,:) ; P(1,:) ; P ; P(npoints,:) ; P(npoints,:) ];
    for i=1:npoints+1
        D3 = (-Q(i,:) + 3*Q(i+1,:) - 3*Q(i+2,:) + Q(i+3,:))/6;
        D2 = ( 3*Q(i,:) - 6*Q(i+1,:) + 3*Q(i+2,:))/6;
        D1 = (-3*Q(i,:) + 3*Q(i+2,:))/6;
        D0 = ( Q(i,:) + 4*Q(i+1,:) + Q(i+2,:))/6;
        t=linspace(0,1,21);
        X1 = D0(1,1)+(D1(1,1)+(D2(1,1)+D3(1,1).*t).*t).*t;
        X2 = D0(1,2)+(D1(1,2)+(D2(1,2)+D3(1,2).*t).*t).*t;
        plot(X1,X2,'b');
        grid on
    end
end
