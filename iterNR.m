function [ M_cur, noconv, NLinc, maxChangeM ] = iterNR( M_cur, M_prev, paramPhys, paramGrid, paramNum, dt, dWmag )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

noconv = 1;
NLinc = 1;
while ( noconv == 1 )
    
    [ F_cur, J_cur ] = calcDiffFieldAtom( M_cur, M_prev, paramPhys, paramGrid, dt, dWmag );
    
    %. numerical derivative
    checkDer = 0;
    if ( checkDer == 1 )
        sdlt = 1e-8;
        J_cur_num = zeros(size(M_cur,1));
        for kk = 1:size(M_cur,1)
            M_cur_c = M_cur;
            M_cur_c(kk) = M_cur(kk) + sdlt;
            F_cur_c = calcDiffFieldAtom( M_cur_c, M_prev, paramPhys, paramGrid, dt, dWmag );
            J_cur_num( :, kk ) = ( F_cur_c - F_cur ) / sdlt;
        end
        verifyJ = max(max(abs(full(J_cur)-J_cur_num)));
    end
    
    Dlt = -J_cur \ F_cur;
    M_cur = M_cur + Dlt;
    
    normF = max(abs(F_cur));
    normDlt = max(abs(Dlt));
    
    if ( normF < paramNum.solveTolEvalF )&&( normDlt < paramNum.solveTolEvalDlt )
        noconv = 0;
    elseif ( NLinc > paramNum.solveMaxIter )
        noconv = 2;
    end
    NLinc = NLinc + 1;
    
end
NLinc = NLinc - 1;

%. maximum change of M locally
changeM = acos( sum( reshape( M_cur.*M_prev, 3, [] ), 1 ) )*180/pi;
maxChangeM = max( changeM );
if ( noconv == 0 )&&( maxChangeM > paramNum.maxAllowedChangeM )
    noconv = 3;
end
   
end

