function [ F_cur, J_cur ] = calcDiffFieldAtom( M_cur, M_prev, paramPhys, paramGrid, dt, dWmag )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

npoints = size(M_cur,1)/3;

linkNumAtoms = paramGrid.linkNumAtoms;
maxLinkNumAtoms = max(linkNumAtoms);

index_i = [ -2 -1 0 -2 -1 0 -2 -1 0 ].';
index_j = [ -2 -2 -2 -1 -1 -1 0 0 0 ].';
J_cur_i = zeros( 9*npoints, maxLinkNumAtoms+1 );
J_cur_j = zeros( 9*npoints, maxLinkNumAtoms+1 );
J_cur_v = zeros( 9*npoints, maxLinkNumAtoms+1 );
F_cur = zeros( 3*npoints, 1 );

alpha = paramPhys.alpha;
beta = paramPhys.beta;
Happ = paramPhys.Happ;
aniKa = paramPhys.aniKa;
aniPa = paramPhys.aniPa;
aniKb = paramPhys.aniKb;
aniPb = paramPhys.aniPb;

for ii = 1:npoints
    
    m_xy_cur = M_cur( (ii*3-2):(ii*3), : );
    m_xy_prev = M_prev( (ii*3-2):(ii*3), : );
    
    %. exchange field
    linkAtoms = paramGrid.linkAtoms_all{ii};
    paramsJ = paramGrid.paramsJ_all{ii};
    indDOF = [ linkAtoms*3-2; linkAtoms*3-1; linkAtoms*3 ];
    M_cur_inter = reshape( M_cur(indDOF(:)), 3, [] );
    M_prev_inter = reshape( M_prev(indDOF(:)), 3, [] );
    Hexc_cur = sum( M_cur_inter.*repmat(paramsJ,3,1), 2 );
    Hexc_prev = sum( M_prev_inter.*repmat(paramsJ,3,1), 2 );
    Hexc = ( Hexc_cur + Hexc_prev )/2;
    
    %. anisotropy field
    Haniso = aniKa * aniPa * ( aniPa.' * ( m_xy_cur + m_xy_prev )/2 ) + aniKb * aniPb * ( aniPb.' * ( m_xy_cur + m_xy_prev )/2 );
    
    %. thermal field
    Hnoise = dWmag( (ii*3-2):(ii*3) )/dt;
    
    %. total
    Heff = Happ + Hexc + Haniso + Hnoise;
    
    %. residual
    s_m_xy_cur = screw(m_xy_cur);
    s_m_xy_prev = screw(m_xy_prev);
    F_cur( (ii*3-2):(ii*3), : ) = m_xy_cur - m_xy_prev - alpha*s_m_xy_prev*m_xy_cur + (beta*dt/2)*(s_m_xy_cur+s_m_xy_prev)*Heff;
    
    %. N.R. derivatives
    d_Hani_d_m = aniKa * ( aniPa * aniPa.' ) / 2 + aniKb * ( aniPb * aniPb.' ) / 2;
    J_cur_add = eye(3) - alpha*s_m_xy_prev - (beta*dt/2)*screw(Heff) + (beta*dt/2)*(s_m_xy_cur+s_m_xy_prev)*d_Hani_d_m;
    J_cur_i( (9*ii-8):(9*ii), 1:(linkNumAtoms(ii)+1) ) = repmat( index_i+ii*3, 1, linkNumAtoms(ii)+1 );
    J_cur_j( (9*ii-8):(9*ii), 1 ) = index_j+ii*3;
    J_cur_j( (9*ii-8):(9*ii), 2:(linkNumAtoms(ii)+1) ) = repmat( index_j, 1, linkNumAtoms(ii) ) + repmat( linkAtoms*3, 9, 1 );
    J_cur_v( (9*ii-8):(9*ii), 1 ) = J_cur_add(:);
    J_cur_add = (beta*dt/4)*(s_m_xy_cur+s_m_xy_prev);
    J_cur_v( (9*ii-8):(9*ii), 2:(linkNumAtoms(ii)+1) ) = repmat( J_cur_add(:), 1, linkNumAtoms(ii) ) .* repmat( paramsJ, 9, 1 );

end

%. assemble Jacobian
rows = J_cur_i(:);
cols = J_cur_j(:);
vals = J_cur_v(:);
zeroElems = (rows==0)|(cols==0);
% J_cur = sparse( rows(~zeroElems), cols(~zeroElems), vals(~zeroElems) );
J_cur = accumarray( [rows(~zeroElems) cols(~zeroElems)], vals(~zeroElems) );

end

