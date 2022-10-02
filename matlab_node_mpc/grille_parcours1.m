function grille_parcours1 = grille_parcours1(xinf,xsup,yinf,ysup,maillage)

grille_parcours1 = [];

for j = 0:fix((xsup-xinf)/maillage)
    
    x = xinf + j * maillage;

    if rem(j,2) == 0

        y = yinf;

        for i = 0:(fix((ysup-yinf)/maillage)-1)
            grille_parcours1= [grille_parcours1 ; x   y + i * maillage];
        end

        grille_parcours1 = [grille_parcours1 ; x   ysup];

    end

    if rem(j,2) == 1

        y = ysup;

        for i = 0:(fix((ysup-yinf)/maillage)-1)
            grille_parcours1  = [grille_parcours1 ; x   y - i * maillage];
        end

        grille_parcours1 = [grille_parcours1 ; x   yinf];

    end

end

 


