
mask1 = (abs(tau_12.*v/norm(m_1-m_2))>1 | abs(tau_34.*v/norm(m_3-m_4))>1);
tau_12(mask1) = [];
tau_34(mask1) = [];
alpha= asin(tau_12.*v/d)*pi;
beta = asin(tau_34.*v/d)*pi;



mask2 = (abs(tau_13.*v/norm(m_1-m_3))>1 | abs(tau_24.*v/norm(m_2-m_4))>1);
tau_13(mask2) = [];
tau_24(mask2) = [];
gamma= asin(tau_13.*v/norm(m_1-m_3))*pi;
delta = asin(tau_24.*v/norm(m_2-m_4))*pi;



mask3 = (abs(tau_23.*v/norm(m_2-m_3))>1 | abs(tau_14.*v/norm(m_1-m_4))>1);
tau_23(mask3) = [];
tau_14(mask3) = [];
epsilon= asin(tau_23.*v/norm(m_2-m_3))*pi;
zeta = asin(tau_14.*v/norm(m_1-m_4))*pi;






tiledlayout(3,3);


function sol = lines_from_points_and_angle(P1, P2, theta)
    validateattributes(P1, {'numeric'}, {'vector','numel',2});
    validateattributes(P2, {'numeric'}, {'vector','numel',2});
    validateattributes(theta, {'numeric'}, {'scalar'});

    x1 = P1(1); y1 = P1(2);
    x2 = P2(1); y2 = P2(2);

    % 1) Line through P1 and P2: a*x + b*y + c = 0
    a1 = y2 - y1;
    b1 = x1 - x2;
    c1 = x2*y1 - x1*y2;

    % 2) Midpoint
    xm = (x1 + x2) / 2;
    ym = (y1 + y2) / 2;

    % 3) Perpendicular bisector (normal = segment direction)
    dx = x2 - x1;
    dy = y2 - y1;
    if dx == 0 && dy == 0
        error('Points P1 and P2 are identical; bisector undefined.');
    end
    a2 = dx;
    b2 = dy;
    c2 = -(a2*xm + b2*ym);

    % 4) Line through midpoint making angle theta with line2
    % direction of line2: v2 = (b2, -a2)
    v2 = [b2, -a2];

    % rotate v2 by theta (matrix multiplication)
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    v3 = (R * v2.').';   % rotated direction as row [vx vy]

    % normal of new line is perpendicular to v3
    a3 = v3(2);
    b3 = -v3(1);
    c3 = -(a3*xm + b3*ym);

    sol = [a1, b1, c1;   % original line
           a2, b2, c2;   % perpendicular bisector
           a3, b3, c3];  % rotated line through midpoint
end









function h = plotlines(D)



xtrue =  linspace(-3,3,2);

[a1, b1, c1] = deal(D(1,1), D(1,2), D(1,3));
[a2, b2, c2] = deal(D(2,1), D(2,2), D(2,3));
[a3, b3, c3] = deal(D(3,1), D(3,2), D(3,3));

% --- Line 1 ---
if abs(b1) > 1e-10
    y_line1 = (-a1*xtrue - c1)/b1;
    %plot(xtrue, y_line1, 'r', 'LineWidth', 2);
else
    % vertical line
    x_line1 = -c1/a1 * ones(size(xtrue));
    %plot(x_line1, xtrue, 'r', 'LineWidth', 2);
end

% --- Line 2 ---
if abs(b2) > 1e-10
    y_line2 = (-a2*xtrue - c2)/b2;
    %plot(xtrue, y_line2, 'g', 'LineWidth', 2);
else
    x_line2 = -c2/a2 * ones(size(xtrue));
    %plot(x_line2, xtrue, 'g', 'LineWidth', 2);
end

% --- Line 3 ---
if abs(b3) > 1e-10
    y_line3 = (-a3*xtrue - c3)/b3;
    h = plot(xtrue, y_line3, 'b', 'LineWidth', 2);
else
    x_line3 = -c3/a3 * ones(size(xtrue));
    h = plot(x_line3, xtrue, 'b', 'LineWidth', 2);
end

end






function [pt, status] = intersect_lines(D1, D2)
% INTERSECT_LINES Intersection of two lines in general form.
% [pt,status] = intersect_lines(D1, D2)
% D1, D2 are 1x3 or 3x1 vectors [a b c] representing a*x + b*y + c = 0
% pt = [x,y] (NaN if no unique intersection)
% status:  1 -> unique intersection
%          0 -> parallel (no intersection)
%         -1 -> coincident (infinite intersections)

    [a1, b1, c1] = deal(D1(3,1), D1(3,2), D1(3,3));
    [a2, b2, c2] = deal(D2(3,1), D2(3,2), D2(3,3));


    A = [a1, b1; a2, b2];
    rhs = -[c1; c2];

    tol = 1e-12;
    detA = a1*b2 - a2*b1;
    if abs(detA) > tol
        xy = A \ rhs;
        pt = xy(:).';
        status = 1;
        return;
    end

    % Lines are (nearly) parallel. Check for coincidence.
    % Two lines a1*x+b1*y+c1=0 and a2*x+b2*y+c2=0 are coincident if
    % the vectors [a b c] are proportional.
    v1 = [a1; b1; c1];
    v2 = [a2; b2; c2];
    if norm(v1) < tol || norm(v2) < tol
        pt = [NaN, NaN];
        status = 0;
        return;
    end

    % Check proportionality by normalizing and comparing
    v1n = v1 / norm(v1);
    v2n = v2 / norm(v2);
    if norm(v1n - sign(v1n.'*v2n)*v2n) < 1e-9   % allow sign flip
        pt = [NaN, NaN];
        status = -1; % coincident
        return;
    end

    % Parallel but not coincident
    pt = [NaN, NaN];
    status = 0;
end





D_11 = zeros(3,3,length(alpha)); 
nexttile;
title("Possible angles with 1-2 pair")
xlim([-3.00 3.00])
ylim([-3.00 3.00])
hold on; 

plot(m_1(1),m_1(2),'ro','MarkerSize',8)
text(m_1(1) - 0.4, m_1(2) + 0.2, 'mic_1')

plot(m_2(1),m_2(2),'ro','MarkerSize',8)
text(m_2(1) + 0.1, m_2(2) + 0.2, 'mic_2')

for k = 1:length(alpha)
    P1 = m_1;
    P2 = m_2;
    theta = alpha;
    D_11(:,:,k) = lines_from_points_and_angle(P1, P2, theta(k));
    line = plotlines(D_11(:,:,k));
end
hold off;






D_21 = zeros(3,3,length(gamma)); 
nexttile;   
xlim([-3.00 3.00])
ylim([-3.00 3.00])
hold on; 
title("Possible angles with 1-3 pair")
plot(m_1(1),m_1(2),'ro','MarkerSize',8)
text(m_1(1) - 0.4, m_1(2) + 0.2, 'mic_1')

plot(m_3(1),m_3(2),'ro','MarkerSize',8)
text(m_3(1) + 0.1, m_3(2) + 0.05, 'mic_3')

plot([m_1(1),m_3(1)], [m_1(2), m_3(2)],'r-' )
for k = 1:length(gamma)
    P1 = m_1;
    P2 = m_3;
    theta = gamma;
    D_21(:,:,k) = lines_from_points_and_angle(P1, P2, theta(k));
    line = plotlines(D_21(:,:,k));
end
hold off;






D_31 = zeros(3,3,length(epsilon));
nexttile;   
xlim([-3.00 3.00])
ylim([-3.00 3.00])
hold on; 
title("Possible angles with 2-3 pair")
plot(m_2(1),m_2(2),'ro','MarkerSize',8)
text(m_2(1) + 0.1, m_2(2) + 0.2, 'mic_2')

plot(m_3(1),m_3(2),'ro','MarkerSize',8)
text(m_3(1) + 0.1, m_3(2) + 0.05, 'mic_3')

plot([m_2(1),m_3(1)], [m_2(2), m_3(2)],'r-' )
for k = 1:length(epsilon)
    P1 = m_2;
    P2 = m_3;
    theta = epsilon;
    D_31(:,:,k) = lines_from_points_and_angle(P1, P2, theta(k));
    line = plotlines(D_31(:,:,k));
end
hold off;






D_12 = zeros(3,3,length(beta)); 
nexttile;   
xlim([-3.00 3.00])
ylim([-3.00 3.00])
hold on; 
title("Possible angles with 3-4 pair")
plot(m_3(1),m_3(2),'ro','MarkerSize',8)
text(m_3(1) + 0.1, m_3(2) + 0.05, 'mic_3')

plot(m_4(1),m_4(2),'ro','MarkerSize',8)
text(m_4(1) + 0.1, m_4(2) - 0.05, 'mic_4')

for k = 1:length(beta)
    P1 = m_3;
    P2 = m_4;
    theta = beta;
    D_12(:,:,k) = lines_from_points_and_angle(P1, P2, theta(k));
    line = plotlines(D_12(:,:,k));
end
hold off;




D_22 = zeros(3,3,length(delta)); 
nexttile;   
xlim([-3.00 3.00])
ylim([-3.00 3.00])
hold on; 
title("Possible angles with 2-4 pair")
plot(m_2(1),m_2(2),'ro','MarkerSize',8)
text(m_2(1) + 0.1, m_2(2) + 0.2, 'mic_2')

plot(m_4(1),m_4(2),'ro','MarkerSize',8)
text(m_4(1) + 0.1, m_4(2) - 0.05, 'mic_4')

plot([m_2(1),m_4(1)], [m_2(2), m_4(2)],'r-' )
for k = 1:length(delta)
    P1 = m_2;
    P2 = m_4;
    theta = delta;
    D_22(:,:,k) = lines_from_points_and_angle(P1, P2, theta(k));
    line = plotlines(D_22(:,:,k));
end
hold off;





D_32 = zeros(3,3,length(zeta)); 
nexttile;   
xlim([-3.00 3.00])
ylim([-3.00 3.00])
hold on; 
title("Possible angles with 1-4 pair")
plot(m_1(1),m_1(2),'ro','MarkerSize',8)
text(m_1(1) - 0.4, m_1(2) + 0.2, 'mic_1')

plot(m_4(1),m_4(2),'ro','MarkerSize',8)
text(m_4(1) + 0.1, m_4(2) - 0.05, 'mic_4')

plot([m_1(1),m_4(1)], [m_1(2), m_4(2)],'r-' )
for k = 1:length(zeta)
    P1 = m_1;
    P2 = m_4;
    theta = zeta;
    D_32(:,:,k) = lines_from_points_and_angle(P1, P2, theta(k));
    line = plotlines(D_32(:,:,k));
end
hold off;




%graph des points




cmap = jet(length(alpha));    % or jet(n), hsv(n), etc. parula

  %clear tout graph present afin de commencer de rien

nexttile
hold on 
title("Drone position with 12-34 pairs")

xlim([-3 3])
ylim([-3 3])
plot(m_1(1),m_1(2),'ro','MarkerSize',8)
text(m_1(1) - 0.4, m_1(2) + 0.2, 'mic_1')
plot(m_2(1),m_2(2),'ro','MarkerSize',8)
text(m_2(1) + 0.1, m_2(2) + 0.2, 'mic_2')
plot(m_3(1),m_3(2),'ro','MarkerSize',8)
text(m_3(1) + 0.1, m_3(2) + 0.05, 'mic_3')
plot(m_4(1),m_4(2),'ro','MarkerSize',8)
text(m_4(1) + 0.1, m_4(2) - 0.05, 'mic_4')


plot([m_1(1),m_2(1)], [m_1(2), m_2(2)],'r-' )
plot([m_3(1),m_4(1)], [m_3(2), m_4(2)],'r-' )

for i = 1:size(D_11,3)
    first = plotlines(D_11(:,:,i));
    second = plotlines(D_12(:,:,i));
    pause(0.25*temps_entre_position)
    [point, status] = intersect_lines(D_11(:,:,i), D_12(:,:,i));
    plot(point(1), point(2), '*', 'Color', cmap(i,:), 'MarkerSize', 8);
    pause(0.15*temps_entre_position)
    delete(first)
    delete(second)
    pause(0.6*temps_entre_position)
end



cmap = jet(length(gamma));    % or jet(n), hsv(n), etc. parula

  

nexttile
hold on 
title("Drone position with 24-13 pairs")

xlim([-3 3])
ylim([-3 3])
plot(m_1(1),m_1(2),'ro','MarkerSize',8)
text(m_1(1) - 0.4, m_1(2) + 0.2, 'mic_1')
plot(m_2(1),m_2(2),'ro','MarkerSize',8)
text(m_2(1) + 0.1, m_2(2) + 0.2, 'mic_2')
plot(m_3(1),m_3(2),'ro','MarkerSize',8)
text(m_3(1) + 0.1, m_3(2) + 0.05, 'mic_3')
plot(m_4(1),m_4(2),'ro','MarkerSize',8)
text(m_4(1) + 0.1, m_4(2) - 0.05, 'mic_4')

plot([m_1(1),m_3(1)], [m_1(2), m_3(2)],'r-' )
plot([m_2(1),m_4(1)], [m_2(2), m_4(2)],'r-' )

for i = 1:size(D_21,3)
    first = plotlines(D_21(:,:,i));
    second = plotlines(D_22(:,:,i));
    pause(0.25*temps_entre_position)
    [point, status] = intersect_lines(D_21(:,:,i), D_22(:,:,i));
    plot(point(1), point(2), '*', 'Color', cmap(i,:), 'MarkerSize', 8);
    pause(0.15*temps_entre_position)
    delete(first)
    delete(second)
    pause(0.6*temps_entre_position)
end





cmap = jet(length(epsilon));    % or jet(n), hsv(n), etc. parula

  

nexttile
hold on
title("Drone position with 14-23 pairs")

xlim([-3 3])
ylim([-3 3])
plot(m_1(1),m_1(2),'ro','MarkerSize',8)
text(m_1(1) - 0.4, m_1(2) + 0.2, 'mic_1')
plot(m_2(1),m_2(2),'ro','MarkerSize',8)
text(m_2(1) + 0.1, m_2(2) + 0.2, 'mic_2')
plot(m_3(1),m_3(2),'ro','MarkerSize',8)
text(m_3(1) + 0.1, m_3(2) + 0.05, 'mic_3')
plot(m_4(1),m_4(2),'ro','MarkerSize',8)
text(m_4(1) + 0.1, m_4(2) - 0.05, 'mic_4')

plot([m_1(1),m_4(1)], [m_1(2), m_4(2)],'r-' )
plot([m_2(1),m_3(1)], [m_2(2), m_3(2)],'r-' )


for i = 1:size(D_31,3)
    first = plotlines(D_31(:,:,i));
    second = plotlines(D_32(:,:,i));
    pause(0.25*temps_entre_position)
    [point, status] = intersect_lines(D_31(:,:,i), D_32(:,:,i));
    plot(point(1), point(2), '*', 'Color', cmap(i,:), 'MarkerSize', 8);
    pause(0.15*temps_entre_position)
    delete(first)
    delete(second)
    pause(0.6*temps_entre_position)
end
