
clear variables;
close all;

fps = 25;
total_t = 6;

%. video output
fn = 'tmp1';
mov = VideoWriter( fn, 'MPEG-4' );
mov.FrameRate = fps;
mov.Quality = 100;
open(mov);

Nframes = round(total_t*fps);

scale = 2.5;
figure( 'Units', 'points', 'Position', [ 50 50 380*scale 230*scale ] );
set( gcf, 'Color', 'w' );

%. data input
load( 'tmp1', 'paramPhys', 'paramGrid', 'coordGrid', 't_all', 'M_all' );
t_min = t_all(1);
t_max = t_all(end);

%. colour scale for m_z
minVal = -0.1;
maxVal = 1;

for ii = 0:Nframes
    
    clf;
    
    hold on;

    t_cur = t_min + (t_max-t_min) * ii/Nframes;
    [ dum, ind ] = min( abs( t_all - t_cur ) );
    
    %. vector plot
    plotField( M_all(:,ind), coordGrid, paramPhys.latSpac, 0.25*scale, minVal, maxVal );
        
    text( 10*scale, 150*scale, [ '{\itt} = ' sprintf( '%.2f', t_cur ) ], 'FontName', 'Segoe UI', 'FontSize', 10*scale, 'Units', 'points', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left' );

    hold off;

    pbaspect( [ 2 1 1 ] );
    limits = [ -0.16 4.16 -0.08 2.08 ];
    axis( limits );
    axis off;

    set( gca, 'Units', 'points' );
    set( gca, 'Position', [ 40*scale 40*scale 300*scale 150*scale ] );
    
    F = getframe(gcf);
    writeVideo( mov, F );
end

close(mov);

