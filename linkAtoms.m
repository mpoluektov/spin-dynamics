function [ paramGrid ] = linkAtoms( coordGrid, paramPhys )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

npoints = size(coordGrid,2);

linkAtoms_all = cell(npoints,1);
paramsJ_all = cell(npoints,1);
Atens_all = zeros(3,npoints);
linkNumAtoms = zeros(1,npoints);

for ii = 1:npoints

    c_xy = coordGrid( :, ii );
    
    %. distances to neighbours
    coordRad = sqrt( ( coordGrid(1,:) - c_xy(1) ).^2 + ( coordGrid(2,:) - c_xy(2) ).^2 );
    coordVec = [ ( coordGrid(1,:) - c_xy(1) ); ( coordGrid(2,:) - c_xy(2) ) ];
    selRad = ( coordRad < paramPhys.cropRad );
    selRad(ii) = 0;
    linkAtoms = find( selRad );
    coordRad_sel = coordRad( linkAtoms );
    coordVec_sel = coordVec( :, linkAtoms );
    
    %. exchange parameters
    paramsJ = paramPhys.coefJ * interLaw( coordRad_sel/paramPhys.latSpac, paramPhys.decCoef );
    
    linkAtoms_all{ii} = linkAtoms;
    paramsJ_all{ii} = paramsJ;
    
    %. coarse-grain exchange
    Ax = ( 1/2 ) * sum( paramsJ .* coordVec_sel(1,:).^2 );
    Ay = ( 1/2 ) * sum( paramsJ .* coordVec_sel(2,:).^2 );
    Axy = ( 1/2 ) * sum( paramsJ .* coordVec_sel(1,:) .* coordVec_sel(2,:) );
    Atens_all(1,ii) = Ax;
    Atens_all(2,ii) = Ay;
    Atens_all(3,ii) = Axy;

    linkNumAtoms(ii) = size(linkAtoms,2);
    
end

paramGrid.linkAtoms_all = linkAtoms_all;
paramGrid.paramsJ_all = paramsJ_all;
paramGrid.Atens_all = Atens_all;
paramGrid.linkNumAtoms = linkNumAtoms;

end

