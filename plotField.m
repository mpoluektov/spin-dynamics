function [ ] = plotField( M_cur, coordGrid, dh, lineWidth, minColVal, maxColVal )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

npoints = size(M_cur,1)/3;

map = [ 0,0,0
        165,0,38
        215,48,39
        244,109,67
        253,174,97
        254,224,139
        255,255,191
        217,239,139
        166,217,106
        102,189,99
        26,152,80
        0,104,55
        150,150,150 ]/255;

for ii = 1:npoints

    %. index -> coordinates
    c_xy = coordGrid( :, ii );

    m_xy = M_cur( (ii*3-2):(ii*3), : );
    
    ang = pi/12;
    relSize = 1/2;
    dc = (maxColVal-minColVal)/11;
    if ( m_xy(3) > maxColVal )
        colStyle = map( 1, : );
    elseif ( m_xy(3) < minColVal )
        colStyle = map( 13, : );
    else
        level = 1+round( 10*(m_xy(3)-minColVal-dc/2)/(maxColVal-minColVal-dc) );
        colStyle = map( 13-level, : );
    end
    plot( [ c_xy(1)-m_xy(1)*dh/2 c_xy(1)+m_xy(1)*dh/2 ], [ c_xy(2)-m_xy(2)*dh/2 c_xy(2)+m_xy(2)*dh/2 ], '-', 'LineWidth', lineWidth, 'Color', colStyle );
    rot1 = [ cos(ang) sin(ang); -sin(ang) cos(ang) ] * m_xy(1:2,:) * dh * relSize;
    rot2 = [ cos(ang) -sin(ang); sin(ang) cos(ang) ] * m_xy(1:2,:) * dh * relSize;
    plot( [ c_xy(1)+m_xy(1)*dh/2 c_xy(1)+m_xy(1)*dh/2-rot1(1) ], [ c_xy(2)+m_xy(2)*dh/2 c_xy(2)+m_xy(2)*dh/2-rot1(2) ], '-', 'LineWidth', lineWidth, 'Color', colStyle );
    plot( [ c_xy(1)+m_xy(1)*dh/2 c_xy(1)+m_xy(1)*dh/2-rot2(1) ], [ c_xy(2)+m_xy(2)*dh/2 c_xy(2)+m_xy(2)*dh/2-rot2(2) ], '-', 'LineWidth', lineWidth, 'Color', colStyle );
    
end

end

