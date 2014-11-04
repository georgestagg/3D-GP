function [mgx,mgy,velx,vely] = gpe2dvelocity(gridx,gridy,ophase)
    
        dspace=(gridx(2)-gridx(1));
        phase = unwrap(ophase);
        dims = size(ophase);
        for i = 2:dims(1)-1
        for j = 2:dims(2)-1
            if (phase(i+1,j)-phase(i-1,j)<-(pi/2.0d0))
                temp1 = phase(i+1,j)-(phase(i-1,j) - pi);
            elseif (phase(i+1,j)-phase(i-1,j)>(pi/2.0d0))
                temp1 = phase(i+1,j)-(phase(i-1,j) + pi);
            else
                temp1 = phase(i+1,j)-phase(i-1,j);
            end
            mgx(i-1,j-1) = gridx(j);
            mgy(i-1,j-1) = gridy(i);
            velx(i-1,j-1) = real(temp1)/dspace;
                    
        end
        end
        phase = unwrap(ophase,[],2);
        for i = 2:dims(1)-1
        for j = 2:dims(2)-1
            if (phase(i,j+1)-phase(i,j-1)<-(pi/2.0d0))
                temp1 = phase(i,j+1)-(phase(i,j-1) - pi);
            elseif (phase(i,j+1)-phase(i,j-1)>(pi/2.0d0))
                temp1 = phase(i,j+1)-(phase(i,j-1) + pi);
            else
                temp1 = phase(i,j+1)-phase(i,j-1);
            end
            vely(i-1,j-1) = real(temp1)/dspace;
        end
        end
        

end
