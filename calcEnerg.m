function [ E_cur ] = calcEnerg( M_cur, paramPhys, paramGrid )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

npoints = size(M_cur,1)/3;

E_cur = zeros( npoints, 1 );

Happ = paramPhys.Happ;
aniKa = paramPhys.aniKa;
aniPa = paramPhys.aniPa;
aniKb = paramPhys.aniKb;
aniPb = paramPhys.aniPb;

for ii = 1:npoints
    
    m_xy_cur = M_cur( (ii*3-2):(ii*3), : );
    
    %. exchange field
    linkAtoms = paramGrid.linkAtoms_all{ii};
    paramsJ = paramGrid.paramsJ_all{ii};
    indDOF = [ linkAtoms*3-2; linkAtoms*3-1; linkAtoms*3 ];
    M_cur_inter = reshape( M_cur(indDOF(:)), 3, [] );
    Hexc_cur = sum( M_cur_inter.*repmat(paramsJ,3,1), 2 );
    paramsJ_sum = sum( paramsJ );
    
    %. energy
    e_exc = -m_xy_cur.' * Hexc_cur/2 + paramsJ_sum/2 ;
    e_aniso = -( aniKa * ( aniPa.' * m_xy_cur )^2 + aniKb * ( aniPb.' * m_xy_cur )^2 )/2 + ( aniKa + aniKb )/2;
    e_app = -m_xy_cur.' * Happ + norm(Happ);
    E_cur(ii) = e_exc + e_aniso + e_app;
    
end

end

