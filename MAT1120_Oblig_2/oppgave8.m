    n = 8;
    t = zeros(n);
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
   C_inv = transpose(C_u);
   S_inv = transpose(S_u);
   
   t_16 = linspace(pi/16, 15*pi/16, 8);
   t = linspace(-pi, pi, 10000);
   f_t16 = (pi.^2. - (t_16.^2)).*exp(t_16./pi);
   f_t16_n = ((pi.^2.- ((-t_16).^2)).*exp((-t_16)./pi));
   f_l = 0.5*(f_t16 + f_t16_n);
   f_o = 0.5*(f_t16 - f_t16_n);
   y_c = mtimes(C_inv, transpose(f_l));
   y_s = mtimes(S_inv, transpose(f_o));
   
   f_16 = y_c(1) + y_c(2)*cos(t) + y_c(3)*cos(2*t) + y_c(4)*cos(3*t) ...
            + y_c(5)*cos(4*t) + y_c(6)*cos(5*t) + y_c(7)*cos(6*t)... 
            + y_c(8)*cos(7*t)+ y_s(1)*sin(t) + y_s(2)*sin(2*t)... 
            + y_s(3)*sin(3*t) + y_s(4)*sin(4*t)+ y_s(5)*sin(5*t)...  
            + y_s(6)*sin(6*t)  + y_s(7)*sin(7*t) + y_s(8)*sin(8*t);
   f = (pi.^2 - t.^2).*exp(t./pi);
   g = figure;
   plot(t, f, ' r ' , t , f_16 , 'b')
   ax = gca;
   ax.FontSize = 15;
   xlabel('x');
   ylabel('y');
   title('Plot av f og f_{16},  f = \pi^2 - t^2 \cdot e^{(t/pi)}');
   %saveas(g,'oppgave_8', 'png')