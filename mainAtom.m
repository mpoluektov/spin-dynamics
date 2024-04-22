function [ ] = mainAtom( filename )
%MAINATOM Atomistic spin dynamics
%   See (Poluektov et al., Comput. Methods Appl. Mech. Eng. 329:219-253, 2018).
%   Code has been written by M. Poluektov.

%% parameters

if ( nargin < 1 )
    error( 'Specify filename' );
end

%. LLG constants (Gilbert form)
paramPhys.alpha = 0;
paramPhys.beta = 1;

%. anisotropy
paramPhys.aniKa = 2*pi^2 + 1;
paramPhys.aniPa = [ 0 1 0 ].';
paramPhys.aniKb = 1;
paramPhys.aniPb = [ 0 0 1 ].';

%. applied field
paramPhys.Happ = [ -1 0 0 ].';

%. interaction parameters
paramPhys.latSpac = 1/8;
paramPhys.cropRad = 8*paramPhys.latSpac + 1e-6;
paramPhys.coefJ = (1/0.926757589165766) / (paramPhys.latSpac^2); 
paramPhys.decCoef = 0.5;

%. initial state parameters
Aexch = 1;
dwInit = 0.85;

%. thermal fluctuations
kTnorm = 0;
paramPhys.coefD = kTnorm * paramPhys.alpha / paramPhys.beta;

%. atomistic grid
Nx = 33;
Ny = 17;

%. times
dt = 0.02;
ninc = 301;

%. tolerances
paramNum.maxAllowedChangeM = 30;
paramNum.solveMaxIter = 20;
paramNum.solveTolEvalF = 1e-8;
paramNum.solveTolEvalDlt = 1e-8;

%% create atomistic grid

%. grid
dx = paramPhys.latSpac;
dy = dx*cos(pi/6);
shiftH = dx/2;
coordGrid = zeros(2,Nx*Ny);
for ii = 1:Ny
    coordGrid(1,(ii*Nx-Nx+1):(ii*Nx)) = (1:Nx)*dx - dx + shiftH*mod(ii-1,2);
    coordGrid(2,(ii*Nx-Nx+1):(ii*Nx)) = (ii-1)*dy*ones(1,Nx);
end

%. void geometry
remCx = 12.5*dx;
remCy = 15.5*dy;
remR = 0.8*paramPhys.latSpac + 1e-6;
coordRem = [ remCx-remR    remCy                ;
             remCx-remR/2  remCy+remR*sqrt(3)/2 ; 
             remCx+remR/2  remCy+remR*sqrt(3)/2 ;
             remCx+remR    remCy                ;
             remCx+remR/2  remCy-remR*sqrt(3)/2 ;
             remCx-remR/2  remCy-remR*sqrt(3)/2 ];

%. remove points - void
npoints = size(coordGrid,2);
indRem = true(1,npoints);
for ii = 1:npoints
    c_xy = coordGrid( :, ii );
    if inpolygon( c_xy(1), c_xy(2), coordRem(:,1), coordRem(:,2) )
        indRem(ii) = false(1);
    end
end
coordGrid = coordGrid(:,indRem);

%% plot atomistic grid

isPlot = 0;
if (isPlot == 1)
    figure(1);
    hold on;
    plot( coordGrid(1,:), coordGrid(2,:), 'ob' );
    plot( [ coordRem(:,1); coordRem(1,1) ], [ coordRem(:,2); coordRem(1,2) ], '-r' );
    axis( [ 0 6 0 6 ] );
    pbaspect( [ 1 1 1 ] );
    hold off;
end

%% create initial state

npoints = size(coordGrid,2);
M_0 = zeros(3*npoints,1);
for ii = 1:npoints
    c_xy = coordGrid( :, ii );

    nu = asin( paramPhys.Happ(1)/paramPhys.aniKa );
    dwPar = sqrt( ( paramPhys.aniKa - paramPhys.aniKb )/Aexch );
    thet = asin(tanh( dwPar*( c_xy(1) - dwInit ) )) + pi/2;
    m_xy = [ sin(nu)           ;
             cos(thet)*cos(nu) ;
             sin(thet)*cos(nu) ];
         
    M_0( (ii*3-2):(ii*3), : ) = m_xy;
end

%% link atoms

paramGrid = linkAtoms( coordGrid, paramPhys );
E_0 = calcEnerg( M_0, paramPhys, paramGrid );

%% initialise

%. variables with calculated data
t_all = (1:ninc)*dt - dt;
M_all = zeros( 3*npoints, ninc );
E_all = zeros( npoints, ninc );
M_all(:,1) = M_0;
E_all(:,1) = E_0;
M_prev = M_0;

%. thermal noise
noise_all = randn( 3*npoints, ninc );

%% calculate 

tic;
for iinc=2:ninc

    %. noise
    dWmag = sqrt(2*paramPhys.coefD) * sqrt(dt) * noise_all( :, iinc ) ;
    
    %. calculate step
    M_cur = M_prev;
    [ M_cur, noconv, NLinc, maxChangeM ] = iterNR( M_cur, M_prev, paramPhys, paramGrid, paramNum, dt, dWmag );
    E_cur = calcEnerg( M_cur, paramPhys, paramGrid );
    
    if ( noconv == 0 )
        %. convergence reached
        M_prev = M_cur;
        M_all(:,iinc) = M_cur;
        E_all(:,iinc) = E_cur;
        
        %. report progress
        if ( mod( iinc, 1 ) == 0 )
            fprintf( '    timestep %.0f/%.0f done, %.0f N.R. increments, maxChangeM = %.2f deg\n', iinc-1, ninc-1, NLinc, maxChangeM );
        end
    else
        %. no convergence
        fprintf( '    no convergence, noconv = %.0f\n', noconv );
        break;
    end
end
elapsed = toc;
fprintf( '    completed, %.0f sec. elapsed\n', elapsed );

%% save

save( filename, 'paramPhys', 'paramGrid', 'paramNum', 'coordGrid', 't_all', 'M_all', 'E_all', 'noise_all' );

end

