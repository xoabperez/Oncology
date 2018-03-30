function [M] = mech_matrix_build(sy,sx,dy,dx,E, nu)
% sy = # of points in the y-direction
% sx = # of points in the x-direction
% dims(1) = grid spacing in the x-direction
% dims(2) = grid spacing in the y-direction
% E = young's modulus
% nu = Poissons ratio

M = zeros(sy*sx*2,sy*sx*2);


    bcf = ones(sy,sx,2);
    count = 0;
    for y = 1:sy
        for x = 1:sx
            count = count + 1;
            itm(y,x) = count;
            if x == 1 || x == sx
                bcf(y,x,1) = 0;
            end
            if y == 1 || y == sy
               bcf(y,x,2) = 0; 
            end               
        end
    end
    
    Ke = E/((1+nu)*(1-2*nu));
    Ks = (1-2*nu)/2;
    
    for y = 1:sy
        for x = 1:sx
            count = itm(y,x);
            county = itm(y,x)+sy*sx;
            
            % X-direction
            if x == 1 || x == sx
                M(count,count) = 1;
            else
                if y == 1
                   % Ux(y,x+1)
                        M(count,itm(y,x+1)) = bcf(y,x+1,1)*(1-nu)*Ke/(dx^2);
                   % Ux(y,x)
                        M(count,itm(y,x)) = -2*Ke*(1-nu)/(dx^2);
                   % Ux(y,x-1)
                        M(count,itm(y,x-1)) = bcf(y,x-1,1)*Ke*(1-nu)/(dx^2);    
                        
                   % Uy(y+1,x+1)
                        M(count,sx*sy+itm(y+1,x+1)) = bcf(y+1,x+1,2)*Ke*(nu+Ks)/(2*dx*dy);
                   % Uy(y,x+1)
                        M(count,sx*sy+itm(y,x+1)) = -bcf(y,x+1,2)*Ke*(nu+Ks)/(2*dx*dy);
                   % Uy(y+!,x-1)
                        M(count,sx*sy+itm(y+1,x-1)) = -bcf(y+1,x-1,2)*Ke*(nu+Ks)/(2*dx*dy);
                   % Uy(y,x-1)
                        M(count,sy*sx+itm(y,x-1)) = bcf(y,x-1,2)*Ke*(nu+Ks)/(2*dx*dy);
                        
                   % Ux(y,x)
                        M(count,itm(y,x)) = M(count,itm(y,x))+bcf(y,x,1)*Ke*(Ks/dy^2);
                   % Ux(y+1,x)
                        M(count,itm(y+1,x)) = -2*bcf(y+1,x,1)*Ke*(Ks/dy^2);
                   % Ux(y+2,x)
                        M(count,itm(y+2,x)) = bcf(y+2,x,1)*Ke*(Ks/dy^2);                
                    
                elseif y == sy
                   % Ux(y,x+1)
                        M(count,itm(y,x+1)) = bcf(y,x+1,1)*(1-nu)*Ke/(dx^2);
                   % Ux(y,x)
                        M(count,itm(y,x)) = -2*Ke*(1-nu)/(dx^2);
                   % Ux(y,x-1)
                        M(count,itm(y,x-1)) = bcf(y,x-1,1)*Ke*(1-nu)/(dx^2);    
                        
                   % Uy(y,x+1)
                        M(count,sx*sy+itm(y,x+1)) = bcf(y,x+1,2)*Ke*(nu+Ks)/(2*dx*dy);
                   % Uy(y-1,x+1)
                        M(count,sx*sy+itm(y-1,x+1)) = -bcf(y-1,x+1,2)*Ke*(nu+Ks)/(2*dx*dy);
                   % Uy(y,x-1)
                        M(count,sx*sy+itm(y,x-1)) = -bcf(y,x-1,2)*Ke*(nu+Ks)/(2*dx*dy);
                   % Uy(y-1,x-1)
                        M(count,sy*sx+itm(y-1,x-1)) = bcf(y-1,x-1,2)*Ke*(nu+Ks)/(2*dx*dy);
                        
                   % Ux(y,x)
                        M(count,itm(y,x)) = M(count,itm(y,x))+bcf(y,x,1)*Ke*(Ks/dy^2);
                   % Ux(y-1,x)
                        M(count,itm(y-1,x)) = -2*bcf(y-1,x,1)*Ke*(Ks/dy^2);
                   % Ux(y-2,x)
                        M(count,itm(y-2,x)) = bcf(y-2,x,1)*Ke*(Ks/dy^2);
                    
                else
                   % Ux(y,x+1)
                        M(count,itm(y,x+1)) = bcf(y,x+1,1)*(1-nu)*Ke/(dx^2);
                   % Ux(y,x)
                        M(count,itm(y,x)) = -2*Ke*(1-nu)/(dx^2);
                   % Ux(y,x-1)
                        M(count,itm(y,x-1)) = bcf(y,x-1,1)*Ke*(1-nu)/(dx^2);    
                        
                   % Uy(y+1,x+1)
                        M(count,sx*sy+itm(y+1,x+1)) = bcf(y+1,x+1,2)*Ke*(nu+Ks)/(4*dx*dy);
                   % Uy(y-1,x+1)
                        M(count,sx*sy+itm(y-1,x+1)) = -bcf(y-1,x+1,2)*Ke*(nu+Ks)/(4*dx*dy);
                   % Uy(y+1,x-1)
                        M(count,sx*sy+itm(y+1,x-1)) = -bcf(y+1,x-1,2)*Ke*(nu+Ks)/(4*dx*dy);
                   % Uy(y-1,x-1)
                        M(count,sy*sx+itm(y-1,x-1)) = bcf(y-1,x-1,2)*Ke*(nu+Ks)/(4*dx*dy);
                        
                   % Ux(y+1,x)
                        M(count,itm(y+1,x)) = bcf(y+1,x,1)*Ke*(1-2*nu)*(1/(2*(dx^2)));
                   % Ux(y,x)
                        M(count,itm(y,x)) = M(count,itm(y,x))-2*bcf(y,x,1)*Ke*(1-2*nu)*(1/(2*(dx^2)));
                   % Ux(y-1,x)
                        M(count,itm(y-1,x)) = bcf(y-1,x,1)*Ke*(1-2*nu)*(1/(2*(dx^2)));
     

                end
   
            end

            
            if y == 1 || y == sy
                M(county,county) = 1;
            else
               if x == 1
                   % Uy(y+1,x)
                        M(county,sx*sy+itm(y+1,x)) = bcf(y+1,x,2)*Ke*(1-nu)/(dy^2);
                   % Uy(y,x)
                        M(county,sx*sy+itm(y,x)) = -2*bcf(y,x,2)*Ke*(1-nu)/(dy^2);                   
                   % Uy(y-1,x)
                        M(county,sx*sy+itm(y-1,x)) = bcf(y-1,x,2)*Ke*(1-nu)/(dy^2);   
                   
                   % Uy(y,x)
                        M(county,sx*sy+itm(y,x)) = M(county,sx*sy+itm(y,x))+1*bcf(y,x,2)*Ke*Ks/(dy^2);                 
                   % Uy(y,x+1)
                        M(county,sx*sy+itm(y,x+1)) = -2*bcf(y,x+1,2)*Ke*Ks/(dy^2); 
                     % Uy(y,x+2)
                        M(county,sx*sy+itm(y,x+2)) = bcf(y,x+2,2)*Ke*Ks/(dy^2);                   
                        
                   %Ux(y+1,x+1)
                        M(county,itm(y+1,x+1)) =  bcf(y+1,x+1,1)*Ke*(nu+Ks)/(2*dx*dy);
                   %Ux(y+1,x-1)
                        M(county,itm(y+1,x)) = -bcf(y+1,x,1)*Ke*(nu+Ks)/(2*dx*dy);               
                   %Ux(y-1,x+1)
                        M(county,itm(y-1,x+1)) = -bcf(y-1,x+1,1)*Ke*(nu+Ks)/(2*dx*dy);                   
                   %Ux(y-1,x-1)
                        M(county,itm(y-1,x)) = bcf(y-1,x,1)*Ke*(nu+Ks)/(2*dx*dy);    
               elseif x == sx
                   % Uy(y+1,x)
                        M(county,sx*sy+itm(y+1,x)) = bcf(y+1,x,2)*Ke*(1-nu)/(dy^2);
                   % Uy(y,x)
                        M(county,sx*sy+itm(y,x)) = -2*bcf(y,x,2)*Ke*(1-nu)/(dy^2);                   
                   % Uy(y-1,x)
                        M(county,sx*sy+itm(y-1,x)) = bcf(y-1,x,2)*Ke*(1-nu)/(dy^2);   
                   
                   % Uy(y,x)
                        M(county,sx*sy+itm(y,x)) = M(county,sx*sy+itm(y,x))+1*bcf(y,x,2)*Ke*Ks/(dy^2);                 
                   % Uy(y,x-1)
                        M(county,sx*sy+itm(y,x-1)) = -2*bcf(y,x-1,2)*Ke*Ks/(dy^2); 
                     % Uy(y,x-2)
                        M(county,sx*sy+itm(y,x-2)) = bcf(y,x-2,2)*Ke*Ks/(dy^2);                   
                        
                   %Ux(y+1,x+1)
                        M(county,itm(y+1,x)) =  bcf(y+1,x,1)*Ke*(nu+Ks)/(2*dx*dy);
                   %Ux(y+1,x-1)
                        M(county,itm(y+1,x-1)) = -bcf(y+1,x-1,1)*Ke*(nu+Ks)/(2*dx*dy);               
                   %Ux(y-1,x+1)
                        M(county,itm(y-1,x)) = -bcf(y-1,x,1)*Ke*(nu+Ks)/(2*dx*dy);                   
                   %Ux(y-1,x-1)
                        M(county,itm(y-1,x-1)) = bcf(y-1,x-1,1)*Ke*(nu+Ks)/(2*dx*dy);    
                   
               else
                   % Uy(y+1,x)
                        M(county,sx*sy+itm(y+1,x)) = bcf(y+1,x,2)*Ke*(1-nu)/(dy^2);
                   % Uy(y,x)
                        M(county,sx*sy+itm(y,x)) = -2*bcf(y,x,2)*Ke*(1-nu)/(dy^2);                   
                   % Uy(y-1,x)
                        M(county,sx*sy+itm(y-1,x)) = bcf(y-1,x,2)*Ke*(1-nu)/(dy^2);   
                   
                   % Uy(y,x+1)
                        M(county,sx*sy+itm(y,x+1)) = bcf(y,x+1,2)*Ke*Ks/(dy^2);
                   % Uy(y,x)
                        M(county,sx*sy+itm(y,x)) = M(county,sx*sy+itm(y,x))-2*bcf(y,x,2)*Ke*Ks/(dy^2);                 
                   % Uy(y,x-1)
                        M(county,sx*sy+itm(y,x-1)) = bcf(y,x-1,2)*Ke*Ks/(dy^2); 
                    
                        
                   %Ux(y+1,x+1)
                        M(county,itm(y+1,x+1)) =  bcf(y+1,x+1,1)*Ke*(nu+Ks)/(4*dx*dy);
                   %Ux(y+1,x-1)
                        M(county,itm(y+1,x-1)) = -bcf(y+1,x-1,1)*Ke*(nu+Ks)/(4*dx*dy);               
                   %Ux(y-1,x+1)
                        M(county,itm(y-1,x+1)) = -bcf(y-1,x+1,1)*Ke*(nu+Ks)/(4*dx*dy);                   
                   %Ux(y-1,x-1)
                        M(county,itm(y-1,x-1)) = bcf(y-1,x-1,1)*Ke*(nu+Ks)/(4*dx*dy);               
                   

               end
                
            end

        end
    end
    M = sparse(M);

end