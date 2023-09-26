function [C_ort1, S_ort2, C_inv, S_inv] = matrix_gen(n)
    t = zeros(n,1);
    C = zeros(n,n);
    S = zeros(n,n);
    for i = 1:n
        t(i) = pi/16 + (i-1)*pi/8;
        for j = 1:n
            C(i,j) = cos((j-1)*t(i));
            S(i,j) = sin(j*t(i));
        end
    end
   C_trans = transpose(C);
   S_trans = transpose(S);
   
   C_ort1 = round(mtimes(C_trans, C));
   S_ort2 = round(mtimes(S_trans, S));  
   
   C_u = zeros(n);
   S_u = zeros(n);
   for k = 1:n
       cj = C(:,k);
       sj = S(:,k);
       C_u(:, k) = cj./(sum(cj.^2));
       S_u(:, k) = sj./(sum(sj.^2));
   end
   C_inv = round(transpose(C_u), 3);
   S_inv = round(transpose(S_u), 3); 
end

