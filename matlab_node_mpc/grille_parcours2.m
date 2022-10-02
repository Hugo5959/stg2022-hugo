function grille_parcours2 = grille_parcours2(xinf,xsup,yinf,ysup,maillage)

grille_parcours2 = [];

for j = 0:(fix((ysup-yinf)/maillage))*2
    
    if rem(j,4) == 0

        x = xinf;
        y = yinf;

        for i = 0:(fix((ysup-yinf)/maillage)-1)
            grille_parcours2= [grille_parcours2 ; x   y + i * maillage];
        end

        if (fix((ysup-yinf)/maillage)) <= 2
            grille_parcours2 = [grille_parcours2 ; x   ysup];
        end

        xinf = xinf + maillage;

    end

    if rem(j,4) == 1

        x = xinf;
        y = ysup;

        for i = 0:(fix((xsup-xinf)/maillage)-1)
            grille_parcours2  = [grille_parcours2 ; x + i * maillage    y];
        end

        if (fix((xsup-xinf)/maillage)) <= 2
            grille_parcours2 = [grille_parcours2 ; xsup   y];
        end

        ysup = ysup - maillage;

    end

    if rem(j,4) == 2

        x = xsup;
        y = ysup;

        for i = 0:(fix((ysup-yinf)/maillage)-1)
            grille_parcours2  = [grille_parcours2 ; x    y - i * maillage];
        end

        if (fix((ysup-yinf)/maillage)) <= 2
            grille_parcours2 = [grille_parcours2 ; x   yinf];
        end

        xsup = xsup - maillage;

    end

    if rem(j,4) == 3

        x = xsup;
        y = yinf;

        for i = 0:(fix((xsup-xinf)/maillage)-1)
            grille_parcours2  = [grille_parcours2 ; x - i * maillage   y];
        end

        if (fix((xsup-xinf)/maillage)) <= 2
            grille_parcours2 = [grille_parcours2 ; xinf   y];
        end

        yinf = yinf + maillage;

    end
    
end