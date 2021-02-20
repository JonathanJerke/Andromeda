/**
 *  coreForce.c
 *
 *
 *  Copyright 2021 Jonathan Jerke and Bill Poirier.
 *  We acknowledge the generous support of Texas Tech University,
 *  the Robert A. Welch Foundation, and the Army Research Office.
 *
 
 *   *   This file is part of Andromeda.
 
 *   *   Andromeda is free software: you can redistribute it and/or modify
 *   *   it under the terms of the GNU General Public License as published by
 *   *   the Free Software Foundation, either version 3 of the License.
 
 *   *   Andromeda is distributed in the hope that it will be useful,
 *   *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   *   GNU General Public License for more details.
 
 *   *   You should have received a copy of the GNU General Public License
 *   *   along with Andromeda.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "coreForce.h"

/**
 * function description
 * @param[out] outString a FILE to print to
 */
void getDescription (   function_label *fn ,double scalar,FILE * outString){
    //param[1] is unused!
    
    if ( fn->fn == nullFunction){
        fprintf(outString,"\tnullFunction\n");
    }else if ( fn->fn == Pseudo ){
        fprintf(outString,"\tPseudo = %1.3f Erf(r/%1.3f)/r \n",scalar*fn->param[0],fn->param[2]);//, fn->param[1]);
    }else if ( fn->fn == Yukawa ){
        fprintf(outString,"\tYukawa = %1.3f exp(- r %1.3f)/r \n",scalar*fn->param[0],fn->param[2]);//, fn->param[1]);
    }else if ( fn->fn == Coulomb ){
        fprintf(outString,"\tCoulomb = %1.3f /r \n",scalar*fn->param[0]);//,fn->param[1]);
    }else if ( fn->fn == Morse ){
        fprintf(outString,"\tMorse = %1.3f (1-exp[-%1.3f *( r - %1.3f )] )^2 -1 \n",scalar*fn->param[0],fn->param[3],fn->param[2]);//,fn->param[1]);
    }else if ( fn->fn == LennardJones ){
        fprintf(outString,"\tLennardJones = %1.3f  rm = %f\n",scalar* fn->param[0],fn->param[2]);
    }
    else if ( fn->fn == Gaussian ){
        fprintf(outString,"\t %f Gaussian --width %f\n",scalar * fn->param[0],fn->param[2]);
    }

    fflush(outString);
}

/**
 * metric description
 * @param[out] outString a string to print
*/
void getMetric (   metric_label mu,FILE * outString){
    switch ( mu.metric ){
        case dirac:
            fprintf(outString,"Dirac @ %f with amp %f\n", mu.beta[0], mu.fn.param[0]);
            break;
        case separateDirac:
            fprintf(outString,"separated Dirac @ %f with amp %f\n", mu.beta[0], mu.fn.param[0]);
            fprintf(outString,"powA %d %d %d\n" , mu.pow[0],mu.pow[1],mu.pow[2]);
            fprintf(outString,"powB %d %d %d\n" , mu.powB[0],mu.powB[1],mu.powB[2]);
            break;
        case interval:
            fprintf(outString,"Interval [%f,%f] with amp %f\n", mu.beta[0],mu.beta[1], mu.fn.param[0]);
        //    printf("deriv %d %d %d\n" , mu.deriv[0],mu.deriv[1],mu.deriv[2]);
        //    printf("pow %d %d %d\n" , mu.pow[0],mu.pow[1],mu.pow[2]);
            break;
        case semiIndefinite:
            fprintf(outString,"Semi-Interval [%f,inf) with amp %f\n", mu.beta[0], mu.fn.param[0]);
        ///    printf("deriv %d %d %d\n" , mu.deriv[0],mu.deriv[1],mu.deriv[2]);
        ///    printf("pow %d %d %d\n" , mu.pow[0],mu.pow[1],mu.pow[2]);
            break;
            default:
                break;

    }
    return;
}

/**
 *Basis overlap Basis
 */
mea BoB (  basisElement_label b1,   basisElement_label b2 ){
    if ( b1.basis == SincBasisElement && b2.basis == SincBasisElement ){
        return SS ( b1.length,b1.length*(b1.index) + b1.origin, b2.length, b2.length*( b2.index ) + b2.origin);
    }else     if ( b1.basis == PeriodicSincBasisElement && b2.basis == PeriodicSincBasisElement ){
        return pSS ( b1.length,b1.length*(b1.index) + b1.origin,b1.grid, b2.length, b2.length*( b2.index ) + b2.origin,b2.grid);
    }

    return 0;
}

/**
 *This differentiates one potential from another
 *
 *The beta dependence is critical
*/
double inverseLaplaceTransform(double beta,   function_label * fl){
    double value2;
    value2 = 0.0;
    if ( fl->fn == Gaussian ){
        return 1;//identity!  differ from Coulomb by some cofactors
    }
    if ( fl->fn == Yukawa ){
        double m = fl->param[2];
        value2  += exp(-(m/2./beta)*(m/2./beta));
    }
    else if ( fl->fn == Morse ){
        double R = fl->param[2];
        double a = fl->param[3];
        value2  += - a * exp (     R * a - (a/beta /2.)*(a/beta /2.))/(beta*beta);
        value2  +=   a * exp ( 2 * R * a - (a/beta    )*(a/beta    ))/(beta*beta);
        value2 *= 1.;
    }
    else if ( fl->fn == Coulomb || fl->fn == Pseudo || fl->fn == nullFunction  ){
            value2 += 1.;//
    } else if ( fl->fn == LennardJones ){
        double rm = fl->param[2];
        value2 +=  1 *   ( 1./60 *  rm * pow(rm * beta,11) ) ;
        value2 += -2 *   ( 1.    *  rm * pow(rm * beta,5 ) ) ;
    }
    return 2./sqrt(pi)*value2;
}

/**
 *Gauss Konrod quadrature
 *
 * https://keisan.casio.com/exec/system/1329114617
 */
double gaussQuad(inta pt , inta nm, inta which ){

    double gk1X [] = {
        0
    };
    double gk1W [] = {
        1
    };

    double gk3X [] = {
        -0.7745966692414833770359,
        0,
        0.7745966692414833770359
    };
    
    double gk3W [] = {
        0.5555555555555555555556,
        0.8888888888888888888889,
        0.555555555555555555556
    };
    
    double gk7X [] = {
        0.949107912342759
        ,0.741531185599394
        ,0.405845151377397
        ,0.0
        ,-0.405845151377397
        ,-0.741531185599394
        ,-0.949107912342759
    };//7
    
    double gk7W [] = {
        0.129484966168870
        ,0.279705391489277
        ,0.381830050505119
        ,0.417959183673469
        ,0.381830050505119
        ,0.279705391489277
        ,0.129484966168870};
    
    
    double gk10X[] = {
        -0.9840853600948424644962,
        -0.9061798459386639927976,
        -0.7541667265708492204408,
        -0.5384693101056830910363,
        -0.2796304131617831934135,
        0,
        0.2796304131617831934135,
        0.5384693101056830910363,
        0.7541667265708492204408,
        0.9061798459386639927976
    };
    
    double gk10W[]= {
        0.04258203675108183286451,
        0.1152333166224733940246,
        0.186800796556492657468,
        0.2410403392286475866999,
        0.272849801912558922341,
        0.2829874178574912132043,
        0.272849801912558922341,
        0.2410403392286475866999,
        0.1868007965564926574678,
        0.1152333166224733940246
    };
    
    
    
    
    double gk15X [] = {
        0.991455371120813
        ,0.949107912342759
        ,0.864864423359769
        ,0.741531185599394
        ,0.586087235467691
        ,0.405845151377397
        ,0.207784955007898
        ,0.0
        ,-0.207784955007898
        ,-0.405845151377397
        ,-0.586087235467691
        ,-0.741531185599394
        ,-0.864864423359769
        ,-0.949107912342759
        ,-0.991455371120813
    };//8
    
    
    double gk15W [] = {
        0.022935322010529
        ,0.063092092629979
        ,0.104790010322250
        ,0.140653259715525
        ,0.169004726639267
        ,0.190350578064785
        ,0.204432940075298
        ,0.209482141084728
        ,0.204432940075298
        ,0.190350578064785
        ,0.169004726639267
        ,0.140653259715525
        ,0.104790010322250
        ,0.063092092629979
        ,0.022935322010529
    };
    
    double gk35X [] = {-0.9984329706060580765167,
        -0.990575475314417335675,
        -0.9746592569674310674486,
        -0.9506755217687677612227,
        -0.9190961368038916732426,
        -0.880239153726985902123,
        -0.834274092850134363076,
        -0.7815140038968014069252,
        -0.722472287372409906023,
        -0.6576711592166907658503,
        -0.5875692123340352515726,
        -0.512690537086476967886,
        -0.4336872952097993710771,
        -0.3512317634538763152972,
        -0.2659465074516820191091,
        -0.1784841814958478558507,
        -0.0895856394252266354625,
        0.00,
        0.0895856394252266354625,
        0.1784841814958478558507,
        0.2659465074516820191091,
        0.3512317634538763152972,
        0.4336872952097993710771,
        0.5126905370864769678863,
        0.5875692123340352515726,
        0.6576711592166907658503,
        0.7224722873724099060232,
        0.7815140038968014069252,
        0.834274092850134363076,
        0.880239153726985902123,
        0.9190961368038916732426,
        0.9506755217687677612227,
        0.9746592569674310674486,
        0.9905754753144173356754,
        0.9984329706060580765167};
    
    double gk35W[] = {0.0042189757937769386559,
        0.0117858375622890860476,
        0.020022233953295123733,
        0.02785672245786342694433,
        0.0352497467518800321001,
        0.04244263020500089117947,
        0.0494361418239590678556,
        0.0559940445300935170357,
        0.062000915268229960367,
        0.067528016718131089923,
        0.072589890114190143846,
        0.0770562230462021636789,
        0.0808366784408178713044,
        0.0839725719912347754276,
        0.0864905322056519368263,
        0.0883087822976044066468,
        0.08936184586778888574733,
        0.0896964219439813653622,
        0.0893618458677888857473,
        0.0883087822976044066468,
        0.086490532205651936826,
        0.08397257199123477542758,
        0.080836678440817871304,
        0.0770562230462021636789,
        0.0725898901141901438458,
        0.067528016718131089923,
        0.062000915268229960367,
        0.05599404453009351703567,
        0.0494361418239590678556,
        0.0424426302050008911795,
        0.0352497467518800321001,
        0.0278567224578634269443,
        0.020022233953295123733,
        0.0117858375622890860476,
        0.00421897579377693865586};
    
    double gk99X [] ={-0.999804199685638727995,
        -0.9988201506066353793618,
        -0.996819814299826469401,
        -0.9937886619441677907601,
        -0.9897642384140644710632,
        -0.984757895914213004359,
        -0.9787573030312066845299,
        -0.97176220090155538014,
        -0.9637900181363554282611,
        -0.954853658674137233555,
        -0.944955059221329473156,
        -0.934100294755810149059,
        -0.922305475453936168755,
        -0.909585655828073285213,
        -0.8959496245835934207352,
        -0.881408445573008910037,
        -0.8659800043151457264443,
        -0.8496821198441657010349,
        -0.832528501487460168517,
        -0.8145344273598554315395,
        -0.7957203207625248381361,
        -0.7761068943454466350181,
        -0.7557118903695143082457,
        -0.734554254237402696214,
        -0.7126570667050088308474,
        -0.6900438244251321135048,
        -0.666735700841571667866,
        -0.642754832419237664057,
        -0.6181268178008792991612,
        -0.5928776941089007124559,
        -0.5670315494953917328631,
        -0.5406132469917260665582,
        -0.5136506284017343823111,
        -0.486171941452492042177,
        -0.4582036915390298548496,
        -0.4297729933415765246586,
        -0.4009095739292798809373,
        -0.3716435012622848888637,
        -0.3420031959918559601078,
        -0.3120175321197487622079,
        -0.2817177082410294178775,
        -0.2511351786125772735072,
        -0.2202997655098053353243,
        -0.1892415924618135864853,
        -0.157992877664368358666,
        -0.126585997269672051068,
        -0.09505164998612337566,
        -0.0634206849826867860288,
        -0.03172586345907315363082,
        0.0,
        0.03172586345907315363082,
        0.06342068498268678602884,
        0.09505164998612337566,
        0.126585997269672051068,
        0.157992877664368358666,
        0.189241592461813586485,
        0.2202997655098053353243,
        0.2511351786125772735072,
        0.2817177082410294178775,
        0.3120175321197487622079,
        0.3420031959918559601078,
        0.3716435012622848888637,
        0.4009095739292798809373,
        0.4297729933415765246586,
        0.4582036915390298548496,
        0.486171941452492042177,
        0.5136506284017343823111,
        0.540613246991726066558,
        0.567031549495391732863,
        0.5928776941089007124559,
        0.618126817800879299161,
        0.6427548324192376640569,
        0.666735700841571667866,
        0.6900438244251321135048,
        0.712657066705008830847,
        0.7345542542374026962137,
        0.7557118903695143082457,
        0.7761068943454466350181,
        0.7957203207625248381361,
        0.81453442735985543154,
        0.8325285014874601685173,
        0.8496821198441657010349,
        0.8659800043151457264443,
        0.881408445573008910037,
        0.8959496245835934207352,
        0.909585655828073285213,
        0.9223054754539361687554,
        0.934100294755810149059,
        0.944955059221329473156,
        0.9548536586741372335552,
        0.9637900181363554282611,
        0.97176220090155538014,
        0.97875730303120668453,
        0.9847578959142130043593,
        0.9897642384140644710632,
        0.9937886619441677907601,
        0.9968198142998264694013,
        0.9988201506066353793618,
        0.999804199685638727995};
    
    double gk99W[] = {0.0005274769683783323143,
        0.0014779639281743620209,
        0.0025218765731496496845,
        0.0035329557014832599892,
        0.0045141331625998620836,
        0.005501420399338078204,
        0.0064998725321664845168,
        0.0074869102764140445198,
        0.0084551856397750456701,
        0.0094175722296862066762,
        0.010378732924116607707,
        0.0113278438578780228795,
        0.0122591748799473589077,
        0.0131792068212079366057,
        0.0140911110472705440377,
        0.0149880989346802956593,
        0.01586572369887289313,
        0.016727898553777318682,
        0.017576872641244826238,
        0.0184077538519825820281,
        0.0192169329826556442117,
        0.0200070643859274929265,
        0.0207798535198561693337,
        0.0215314818077816867538,
        0.022258913930129946636,
        0.0229641216311680166298,
        0.0236484943738631350276,
        0.02430890292981194970869,
        0.0249427312116380137806,
        0.0255515676298983967962,
        0.0261366295727170561997,
        0.0266952731071580377427,
        0.027225205808796711996,
        0.0277278075508688728381,
        0.0282042205700292005214,
        0.0286521667109019998496,
        0.0290696140121504756931,
        0.0294578454056420734076,
        0.0298179973237107521896,
        0.030148081066799933145,
        0.0304462797396858877386,
        0.030713855857739083392,
        0.030951992615566392528,
        0.03115893749659741587,
        0.0313330525912310597589,
        0.03147563582699890351852,
        0.031587959458685973161,
        0.03166846605661171584462,
        0.031715664503973344041,
        0.0317309313985218944091,
        0.031715664503973344041,
        0.0316684660566117158446,
        0.031587959458685973161,
        0.031475635826998903519,
        0.0313330525912310597589,
        0.03115893749659741587,
        0.03095199261556639252752,
        0.030713855857739083392,
        0.0304462797396858877386,
        0.030148081066799933145,
        0.0298179973237107521896,
        0.029457845405642073408,
        0.0290696140121504756931,
        0.0286521667109019998496,
        0.02820422057002920052136,
        0.027727807550868872838,
        0.027225205808796711996,
        0.026695273107158037743,
        0.0261366295727170561997,
        0.025551567629898396796,
        0.0249427312116380137806,
        0.0243089029298119497087,
        0.0236484943738631350276,
        0.02296412163116801663,
        0.02225891393012994663599,
        0.0215314818077816867538,
        0.020779853519856169334,
        0.0200070643859274929265,
        0.0192169329826556442117,
        0.0184077538519825820281,
        0.017576872641244826238,
        0.016727898553777318682,
        0.01586572369887289313,
        0.014988098934680295659,
        0.0140911110472705440377,
        0.0131792068212079366057,
        0.012259174879947358908,
        0.0113278438578780228795,
        0.010378732924116607707,
        0.0094175722296862066762,
        0.0084551856397750456701,
        0.0074869102764140445198,
        0.0064998725321664845168,
        0.005501420399338078204,
        0.0045141331625998620836,
        0.0035329557014832599892,
        0.002521876573149649685,
        0.0014779639281743620209,
        0.00052747696837833231426};
    double *gkX,*gkW;
    inta ngk;
    if ( pt ==1 ){
        gkX = gk1X;
        gkW = gk1W;
        ngk = 1;
    }
    if ( pt == 3 ){
        gkX = gk3X;
        gkW = gk3W;
        ngk = 3;
    }else
    if ( pt == 7 ) {
        gkX = gk7X;
        gkW = gk7W;
        ngk = 7;
    }
        else if ( pt == 10 ){
            gkX = gk10X;
            gkW = gk10W;
            ngk = 10;
        
    }else if ( pt == 15 ){
        gkX = gk15X;
        gkW = gk15W;
        ngk = 15;
    }else if ( pt == 35 ){
        gkX = gk35X;
        gkW = gk35W;
        ngk = 35;
    }else if ( pt == 99 ){
        gkX = gk99X;
        gkW = gk99W;
        ngk = 99;
    }else {
        exit(0);
    }
    
    if ( which )
        return 0.5*(1+gkX[ngk-1-nm]);
    else
        return 0.5*gkW[ngk-1-nm];//shift to interval [0,1]
}


/**
 *building quantum operators for oneBody and twoBody interactions
 *
 */
inta separateInteraction(   sinc_label *f,double scalar, double * position,inta invert,inta act,  blockType bl, double adjustOne,    division load,  metric_label metric,  spinType spin,  division basis ,inta particle1,  bodyType body,inta embed){
      genusType hidden;
      sinc_label f1 = *f;
      division temp , currLoop, currChain,newLabel;
    
    double oneL,twoL;
    inta perm[7],op[7];
    
    inta space1 = 0;
    {
        inta space;
        for ( space = 0 ;space < SPACE  ; space++)
            if ( f1.canon[space].body != nada )
                if ( particle1 == f1.canon[space].label ){
                    space1 = space;
                    break;
                }
        
    }
    
    {
          division li = load;
        while ( f1.name[li].chainNext != nullName)
            li =f1.name[li].chainNext;
        currChain = li;
        temp = eikonBuffer;
    }
    if ( metric.fn.fn == nullFunction){
        printf("null func\n");
        return 0;
    }
    inta spacy, n1[SPACE];
    length1(f1,n1);
    
    inta i,beta,I1,space,I2,N1;
    
    double constant,x,g;
        
    inta section=2,si,ngk, intv = metric.fn.interval,flagConstants=0;
    
    if ( metric.metric == interval )
        section = 0;
    if ( metric.metric == semiIndefinite)
        section = 1;
    if ( metric.metric == dirac)
        section = 2;
    
    if ( metric.metric == pureSemiIndefinite){
        section = 1;
        flagConstants = 1;
    }
    if ( metric.metric == pureInterval){
        flagConstants = 1;
        section = 0;
    }
    if ( section < 2 ){
        ngk = intv;
    }
    else {
        ngk = 1;
    }

    for ( beta = 0; beta < ngk ; beta++){//beta is an index.
        if ( section == 1 ){
            g = gaussQuad(ngk,beta,1);// [1, inf)
            constant = gaussQuad(ngk, beta, 0);
            
            x = ( g ) / (1. - g)+1 ;
            constant /= (1.-g)*(1.-g);
            
            
            x *=  metric.beta[0];
            constant *= metric.beta[0];
            if ( ! flagConstants )
                constant *= inverseLaplaceTransform(x,&metric.fn)*scalar;
            else
                constant = 1;
        }else if ( section == 0 ) {
            g = gaussQuad(ngk,beta,1);//interval [0,1]
            constant = gaussQuad(ngk, beta, 0);
            
            x = g;
            
            
            x *=  (metric.beta[1]-metric.beta[0]);
            constant *= (metric.beta[1]-metric.beta[0]);
            x += metric.beta[0];
            if ( ! flagConstants )
                constant *= inverseLaplaceTransform(x,&metric.fn)*scalar;
            else
                constant = 1;

        } else {
            x = metric.beta[0];// value;
            constant = scalar;
        }
        //printf("x %f \n const %f\n",x,constant);
        tClear(f1,temp);
        tId(f1,temp,0);
        
        //x is beta.
        tClear(f1,temp);
        zero(f1,temp,0);
        inta invertSign;
    
        newLabel = anotherLabel(f,0,nada);
        f1.name[currChain].chainNext = newLabel;
        f1.name[newLabel].species = eikon;
        f1.name[newLabel].multId = beta;
        f1.name[newLabel].spinor = spin;

        ///all equal-beta chained Ops will multiply on each beta index. i.e. H2+
        currChain = newLabel;
        currLoop = currChain;
        
        
        for ( hidden = eikonDiagonal ; hidden <= eikonDiagonal + imin(body,metric.fn.contr);hidden++ )
            {
                
                double oneOri,twoOri,grpL;
                invertSign = 1;

            for ( space = 0 ;space < SPACE  ; space++)
                if ( f1.canon[space].body != nada )
                    if ( f1.canon[space].label == particle1 )
                {
                    if ( body == one ){
                        commandSA(f1.canon[space].body, f1.name[newLabel].space[space].act,tv1 , bl, perm, op);

                            oneL = f1.canon[space].particle[op[0]+1].lattice;
                            oneOri = f1.canon[space].particle[op[0]+1].origin;
                            ///position of left edge...

                            double * te = streams(f1, temp, 0, space);
                            N1 = n1[space];

                                for ( si = 0 ; si < N1; si++){
                                        I1 = si;//
                                    te[si] = momentumIntegralInTrain(x*oneL, ((I1*oneL+oneOri)-position[f1.canon[space].space])/oneL,1, hidden, body);
                                    if ( invertSign  ){
                                            te[si] *= constant;
                                            for ( spacy = 0 ; spacy < embed ; spacy++)
                                                te[si] *= momentumIntegralInTrain(x*oneL, 0,1, hidden, body);
                                    }
                                    if ( alloc(f1, temp, space) < si ){
                                        printf("creation of oneBody, somehow allocations of vectors are too small. %d\n",newLabel);
                                        exit(0);
                                    }

                                }

                        
                        
                        
                    }else
                    ///conditional body 2
                        if  ( body == two )
                    {
                        commandSA(f1.canon[space].body, f1.name[newLabel].space[space].act,e12 , bl, perm, op);
                        oneL = f1.canon[space].particle[op[0]+1].lattice*adjustOne;
                        oneOri = f1.canon[space].particle[op[0]+1].origin*adjustOne;
                        twoL = f1.canon[space].particle[op[1]+1].lattice;
                        twoOri = f1.canon[space].particle[op[1]+1].origin;
                        
                        if ( oneL > twoL ){
                            grpL = fabs(oneL);
                        }else {
                            grpL = fabs(twoL);

                        }
                        N1 = n1[space];

                        double * te = streams(f1, temp, 0, space);

                                si = 0;
                                for ( I2 = 0; I2 < N1; I2++)
                                    for ( I1 = 0 ; I1 < N1; I1++)
                                         {
                                             te[si] = momentumIntegralInTrain(x*grpL, ((oneL*I1+oneOri)-(twoL*I2+twoOri))/grpL,1, hidden, body);
                                         if ( invertSign ){
                                             te[si] *= constant;
                                             for ( spacy = 0 ; spacy < embed ; spacy++)
                                                 te[si] *= momentumIntegralInTrain(x*grpL, 0,1, hidden, body);
                                           }
                                         }
                                         if ( isnan(te[si])){
                                             printf("error in matrix coreForce.c");
                                             exit(0);
                                         }

                                             
                                        si++;
                                         if ( alloc(f1, temp, space) < si ){
                                             printf("creation of twoBody, somehow allocations of vectors are too small. %d\n",newLabel);
                                             exit(0);
                                         }
                            }
                    invertSign = 0;
                }
                
            newLabel = anotherLabel(f,particle1,body);
            ///if this is the first of two entries
            f1.name[newLabel].multId = 0;
            f1.name[newLabel].spinor = spin;
            for ( space = 0 ;space < SPACE  ; space++)
                if ( f1.canon[space].body != nada ){
                    f1.name[newLabel].space[space].act = act;
                    if ( f1.canon[space].label == particle1 )
                            f1.name[newLabel].space[space].body = body;
                    
            }
            
            tEqua(f1, newLabel, 0, temp, 0);
            
            f1.name[currLoop].loopNext = newLabel;
            f1.name[newLabel].species = hidden;
            currLoop = newLabel;
            }
    }
    return 0;
}

/**
 *building quantum operators for oneBody and twoBody interactions
 *Novel inputs:
 *--shape of external field,i.e., flat , Sine , or Gaussian
 *--Novel metric on momentum integral as sufficient
 *---basically a split operator has all kinds of degrees of internal freedom, specifiy
 *
 */
inta periodicInteraction( sinc_label *f,double scalar, double * position,inta invert,inta act,  blockType bl, division load,  metric_label metric,  spinType spin, momentumIntegralSpecs *specs, division basis ,inta particle1,  bodyType body){
    
    
    sinc_label f1 = *f;
      division temp , currLoop, currChain,newLabel;
    
    double oneL,twoL;
    inta perm[7],op[7];
    
    
    
    {
          division li = load;
        while ( f1.name[li].chainNext != nullName)
            li =f1.name[li].chainNext;
        currChain = li;
        temp = eikonBuffer;
    }
    if ( metric.fn.fn == nullFunction){
        printf("null func\n");
        return 0;
    }
    inta  n1[SPACE];
    length1(f1,n1);
    
    inta beta,I1,space,I2,N1;
    
    double constant,x,g;
        
    inta section=2,si,ngk, intv = metric.fn.interval,flagConstants=0;
    
    if ( metric.metric == interval )
        section = 0;
    if ( metric.metric == semiIndefinite)
        section = 1;
    if ( metric.metric == dirac)
        section = 2;
    
    if ( metric.metric == pureSemiIndefinite){
        section = 1;
        flagConstants = 1;
    }
    if ( metric.metric == pureInterval){
        flagConstants = 1;
        section = 0;
    }
    if ( section < 2 ){
        ngk = intv;
    }
    else {
        ngk = 1;
    }

    for ( beta = 0; beta < ngk ; beta++){//beta is an index.
        if ( section == 1 ){
            g = gaussQuad(ngk,beta,1);// [1, inf)
            constant = gaussQuad(ngk, beta, 0);
            
            x = ( g ) / (1. - g)+1 ;
            constant /= (1.-g)*(1.-g);
            
            
            x *=  metric.beta[0];
            constant *= metric.beta[0];
            if ( ! flagConstants )
                constant *= inverseLaplaceTransform(x,&metric.fn)*scalar;
            else
                constant = 1;
        }else if ( section == 0 ) {
            g = gaussQuad(ngk,beta,1);//interval [0,1]
            constant = gaussQuad(ngk, beta, 0);
            
            x = g;
            
            
            x *=  (metric.beta[1]-metric.beta[0]);
            constant *= (metric.beta[1]-metric.beta[0]);
            x += metric.beta[0];
            if ( ! flagConstants )
                constant *= inverseLaplaceTransform(x,&metric.fn)*scalar;
            else
                constant = 1;

        } else {
            x = metric.beta[0];// value;
            constant = scalar;
        }
        //x is beta.
        tClear(f1,temp);
        zero(f1,temp,0);
        inta invertSign;
    
        newLabel = anotherLabel(f,0,nada);
        f1.name[currChain].chainNext = newLabel;
        f1.name[newLabel].species = eikon;
        f1.name[newLabel].multId = beta;
        f1.name[newLabel].spinor = spin;

        ///all equal-beta chained Ops will multiply on each beta index. i.e. H2+
        currChain = newLabel;
        currLoop = currChain;
        
            
        inta momentumIndex,momentumLength, bodyIndex ;
                
        
        for ( bodyIndex = 0 ; bodyIndex < body ; bodyIndex++){
            if ( specs[bodyIndex].metric == zeroMomentum )
                momentumLength = 0;
            else if ( specs[bodyIndex].metric == discreteMomentum )
                momentumLength = specs[bodyIndex].interval;
    
            for ( momentumIndex = -momentumLength ;  momentumIndex <= momentumLength ; momentumIndex++)
            {

                
                double iL,oneOri,twoOri,iO;
                invertSign = 1;

            for ( space = 0 ;space < SPACE  ; space++)
                if ( f1.canon[space].body != nada )
                    if ( f1.canon[space].label == particle1 )
                    {

                        commandSA(f1.canon[space].body, f1.name[newLabel].space[space].act,e12 , bl, perm, op);
                        oneL = f1.canon[space].particle[op[0]+1].lattice;
                        oneOri = f1.canon[space].particle[op[0]+1].origin;
                            
                        if ( body == two) {
                            twoL = f1.canon[space].particle[op[1]+1].lattice;
                            twoOri = f1.canon[space].particle[op[1]+1].origin;
                        }else {
                            twoL= 0.;
                            twoOri = 0.;
                        }
                        N1 = n1[space];


                        DCOMPLEX tc;
                        double * te = streams(f1, temp, 0, space);
                        double * tec = streams(f1, temp, 1, space);

                        si = 0;
                        for ( I2 = 0; I2 < N1; I2++)
                            for ( I1 = 0 ; I1 < N1; I1++)
                                     {
                                             
                                             if ( bodyIndex == 0){
                                                 iL = oneL;
                                                 iO = oneOri;
                                             }
                                             else{
                                                 iL = twoL;
                                                 iO = twoOri;
                                             }
                                         if (specs[bodyIndex].opQ)
                                             tc = periodicSincfourierIntegralInTrain( (I1*iL+iO)/iL/(N1),
                                                                                  (I2*iL+iO)/iL/(N1),  N1,  (1-2*bodyIndex)*momentumIndex);
                                         else
                                             tc = 1.0;
                                         
                                            if ( invertSign && bodyIndex == 0 ){
                                                double gaussianKernel = exp(-pow(pi*momentumIndex/(x*iL),2.))/2./sqrt(pi)/(x*iL);
                                                ///multiply of Gaussian here one first particle.
                                                tc *= constant*gaussianKernel;
                                                
                                                ///for periodic-dirac located external fields
                                                if ( body == one && specs[1].opQ == 0 )
                                                    tc *= cexp(-I*position[f1.canon[space].space]*momentumIndex/iL);
                                            }
                                         te[si] = creal(tc);
                                         tec[si] = cimag(tc);

                                         
                                            if ( isnan(creal(tc)) || isnan(cimag(tc))){
                                                printf("periodicInteraction error\n");
                                                exit(0);
                                            }
                                         si++;
                                         if ( alloc(f1, temp, space) < si ){
                                             printf("creation of twoBody, somehow allocations of vectors are too small. %d\n",newLabel);
                                             exit(0);
                                         
                                     }
                                }
                            }
                    invertSign = 0;
                
                
            newLabel = anotherLabel(f,particle1,body);
            ///if this is the first of two entries
            f1.name[newLabel].multId = 0;
            if ( bodyIndex == 0 && body == two)
                f1.name[newLabel].multId = 1;
            f1.name[newLabel].spinor = spin;
            for ( space = 0 ;space < SPACE  ; space++)
                if ( f1.canon[space].body != nada ){
                    f1.name[newLabel].space[space].act = act;
                    if ( f1.canon[space].label == particle1 ){
                        if ( f1.canon[space].basis == SincBasisElement )
                            f1.name[newLabel].space[space].body = body;
                        else
                            f1.name[newLabel].space[space].body = one;
                        if ( f1.canon[space].basis == SincBasisElement )
                            f1.name[newLabel].space[space].block = bl;
                        else if ( f1.canon[space].basis == PeriodicSincBasisElement ){
                            switch( bl ){
                                case e12:
                                    if ( bodyIndex == 0 )
                                        f1.name[newLabel].space[space].block= tv1;
                                    else
                                        f1.name[newLabel].space[space].block= tv2;
                                break;
                                case e13:
                                    if ( bodyIndex == 0 )
                                        f1.name[newLabel].space[space].block= tv1;
                                    else
                                        f1.name[newLabel].space[space].block= tv3;
                                    break;
                                case e23:
                                    if ( bodyIndex == 0 )
                                        f1.name[newLabel].space[space].block= tv2;
                                    else
                                        f1.name[newLabel].space[space].block= tv3;
                                    break;
                                case e14:
                                    if ( bodyIndex == 0 )
                                        f1.name[newLabel].space[space].block= tv1;
                                    else
                                        f1.name[newLabel].space[space].block= tv4;
                                    break;
                                case e24:
                                    if ( bodyIndex == 0 )
                                        f1.name[newLabel].space[space].block= tv2;
                                    else
                                        f1.name[newLabel].space[space].block= tv4;
                                    break;
                                case e34:
                                    if ( bodyIndex == 0 )
                                        f1.name[newLabel].space[space].block= tv3;
                                    else
                                        f1.name[newLabel].space[space].block= tv4;
                                    break;
                                default:
                                    f1.name[newLabel].space[space].block= bl;
                                ///distribute bl commands
                            }
                        }
                    }
            }
            {
                inta n;
                for ( n = 0 ; n < spin ; n++)
                    tEqua(f1, newLabel, n, temp, n);
            }
            f1.name[currLoop].loopNext = newLabel;
            f1.name[newLabel].species = eikonSplit;
            currLoop = newLabel;
            }
        }
    }
    return 0;
}


/**
 *building quantum operators for oneBody and twoBody interactions
 *For GaussianSinc Basis
*/
inta factoredInteraction(   sinc_label *f,double scalar, double * position,inta invert, blockType bl, double adjustOne,    division load,  metric_label metric,  spinType spin,  division basis ,inta particle1,  bodyType body,inta embed){
      genusType hidden;
      sinc_label f1 = *f;
      division temp , currLoop, currChain,newLabel;
    
    double oneL,twoL;
    inta perm[7],op[7];
    
    inta space1 = 0;
    {
        inta space;
        for ( space = 0 ;space < SPACE  ; space++)
            if ( f1.canon[space].body != nada )
                if ( particle1 == f1.canon[space].label ){
                    space1 = space;
                    break;
                }
        
    }
    
    {
          division li = load;
        while ( f1.name[li].chainNext != nullName)
            li =f1.name[li].chainNext;
        currChain = li;
        temp = eikonBuffer;
    } 
    if ( metric.fn.fn == nullFunction){
        printf("null func\n");
        return 0;
    }
    inta n1[SPACE];
    length1(f1,n1);
    
    inta beta,I1,space,N1;
    
    double constant,x,g;
        
    inta section=2,si,ngk, intv = metric.fn.interval,flagConstants=0;
    
    if ( metric.metric == interval )
        section = 0;
    if ( metric.metric == semiIndefinite)
        section = 1;
    if ( metric.metric == dirac)
        section = 2;
    
    if ( metric.metric == pureSemiIndefinite){
        section = 1;
        flagConstants = 1;
    }
    if ( metric.metric == pureInterval){
        flagConstants = 1;
        section = 0;
    }
    if ( section < 2 ){
        ngk = intv;
    }
    else {
        ngk = 1;
    }

    for ( beta = 0; beta < ngk ; beta++){//beta is an index.
        if ( section == 1 ){
            g = gaussQuad(ngk,beta,1);// [1, inf)
            constant = gaussQuad(ngk, beta, 0);
            
            x = ( g ) / (1. - g)+1 ;
            constant /= (1.-g)*(1.-g);
            
            
            x *=  metric.beta[0];
            constant *= metric.beta[0];
            if ( ! flagConstants )
                constant *= inverseLaplaceTransform(x,&metric.fn)*scalar;
            else
                constant = 1;
        }else if ( section == 0 ) {
            g = gaussQuad(ngk,beta,1);//interval [0,1]
            constant = gaussQuad(ngk, beta, 0);
            
            x = g;
            
            
            x *=  (metric.beta[1]-metric.beta[0]);
            constant *= (metric.beta[1]-metric.beta[0]);
            x += metric.beta[0];
            if ( ! flagConstants )
                constant *= inverseLaplaceTransform(x,&metric.fn)*scalar;
            else
                constant = 1;

        } else {
            x = metric.beta[0];// value;
            constant = scalar;
        }
        
        //x is beta.
        tClear(f1,temp);
        zero(f1,temp,0);
        inta invertSign;
    
        newLabel = anotherLabel(f,0,nada);
        f1.name[currChain].chainNext = newLabel;
        f1.name[newLabel].species = eikon;
        f1.name[newLabel].multId = beta;
        f1.name[newLabel].spinor = spin;

        ///all equal-beta chained Ops will multiply on each beta index. i.e. H2+
        currChain = newLabel;
        currLoop = currChain;
        
        
        for ( hidden = eikonDiagonal ; hidden <= eikonDiagonal + imin(body,metric.fn.contr);hidden++ )
            {
                inta momentumIndex,momentumLength, bodyIndex ;
                momentumLength = 2*f1.canon[space1].count1Basis;
                for ( momentumIndex = -momentumLength ;  momentumIndex <= momentumLength ; momentumIndex++)
                    for ( bodyIndex = 0 ; bodyIndex < body ; bodyIndex++)
                {
                
                
                    double iL,oneOri,twoOri,grpL,iO;
                    invertSign = 1;

                for ( space = 0 ;space < SPACE  ; space++)
                    if ( f1.canon[space].body != nada )
                        if ( f1.canon[space].label == particle1 )
                    {
                    ///conditional body 2
                        {
                            commandSA(f1.canon[space].body, f1.name[newLabel].space[space].act,e12 , bl, perm, op);
                            oneL = f1.canon[space].particle[op[0]+1].lattice*adjustOne;
                            oneOri = f1.canon[space].particle[op[0]+1].origin*adjustOne;
                            
                            if ( body == two) {
                                twoL = f1.canon[space].particle[op[1]+1].lattice;
                                twoOri = f1.canon[space].particle[op[1]+1].origin;
                            }else {
                                twoL= 0.;
                                twoOri = 0.;
                            }
                            if ( oneL > twoL ){
                                grpL = fabs(oneL);
                            }else {
                                grpL = fabs(twoL);

                            }
                            N1 = n1[space];


                            DCOMPLEX tc;
                            double * te = streams(f1, temp, 0, space);
                            double * tec = streams(f1, temp, 1, space);

                            si = 0;
                            for ( I1 = 0 ; I1 < N1; I1++)
                                {
                                             
                                             ///HERE,  assign a quadature pattern
                                                 double momentum = (momentumIndex);
                                                 if ( bodyIndex == 0){
                                                     iL = oneL;
                                                     iO = oneOri;
                                                 }
                                                 else{
                                                     iL = twoL;
                                                     iO = twoOri;
                                                 }
                                                 ///HERE somehow need to populate gamma-s
                                                 ///ALSO NEED a Metric
                                                 
                                                 //tc = gaussianSincfourierIntegralInTrain(iL, <#double gammax#>, <#double x#>, <#double gammay#>, <#double y#>, si, momentum, hidden);
                                                 te[si] = creal(tc);
                                                 tec[si] = cimag(tc);
                                                 
                                                if ( invertSign && bodyIndex == 0 ){
                                                    
                                                    ///find a quadrature given
                                                    ///momentumIndex
                                                    ///
                                                    double gaussianKernel = 1;
                                                    
                                                    ///multiply of Gaussian here one first particle.
                                                    te[si] *= constant*gaussianKernel;
                                                    tec[si] *= constant*gaussianKernel;
                                                }
                                             
                                             if ( isnan(te[si]))
                                                 printf("%f\n",te[si]);
                                             si++;
                                             if ( alloc(f1, temp, space) < si ){
                                                 printf("creation of twoBody, somehow allocations of vectors are too small. %d\n",newLabel);
                                                 exit(0);
                                             }
                                         
                                    }
                                }
                    invertSign = 0;
                }
                
            newLabel = anotherLabel(f,particle1,body);
            ///if this is the first of two entries
            f1.name[newLabel].multId = 0;
            if ( f1.canon[space].basis == PeriodicSincBasisElement || f1.canon[space].basis == GaussianSincBasisElement )
                if ( bodyIndex == 0 && body == two)
                    f1.name[newLabel].multId = 1;
            f1.name[newLabel].spinor = spin;
            for ( space = 0 ;space < SPACE  ; space++)
                if ( f1.canon[space].body != nada ){
                    f1.name[newLabel].space[space].act = 1;
                    if ( f1.canon[space].label == particle1 ){
                        if ( f1.canon[space].basis == SincBasisElement )
                            f1.name[newLabel].space[space].body = body;
                        else
                            f1.name[newLabel].space[space].body = one;
                        if ( f1.canon[space].basis == SincBasisElement )
                            f1.name[newLabel].space[space].block = bl;
                        else if ( f1.canon[space].basis == PeriodicSincBasisElement ){
                            switch( bl ){
                                case e12:
                                    if ( bodyIndex == 0 )
                                        f1.name[newLabel].space[space].block= tv1;
                                    else
                                        f1.name[newLabel].space[space].block= tv2;
                                break;
                                case e13:
                                    if ( bodyIndex == 0 )
                                        f1.name[newLabel].space[space].block= tv1;
                                    else
                                        f1.name[newLabel].space[space].block= tv3;
                                    break;
                                case e23:
                                    if ( bodyIndex == 0 )
                                        f1.name[newLabel].space[space].block= tv2;
                                    else
                                        f1.name[newLabel].space[space].block= tv3;
                                    break;
                                case e14:
                                    if ( bodyIndex == 0 )
                                        f1.name[newLabel].space[space].block= tv1;
                                    else
                                        f1.name[newLabel].space[space].block= tv4;
                                    break;
                                case e24:
                                    if ( bodyIndex == 0 )
                                        f1.name[newLabel].space[space].block= tv2;
                                    else
                                        f1.name[newLabel].space[space].block= tv4;
                                    break;
                                case e34:
                                    if ( bodyIndex == 0 )
                                        f1.name[newLabel].space[space].block= tv3;
                                    else
                                        f1.name[newLabel].space[space].block= tv4;
                                    break;
                                default:
                                    f1.name[newLabel].space[space].block= bl;
                                ///distribute bl commands
                            }
                        }
                    }
            }
                {
                            inta n;
                            for ( n = 0 ; n < spin ; n++){
                                f1.name[temp].Current[n]= 1;
                                tEqua(f1, newLabel, n, temp, n);
                            }
                }
                f1.name[currLoop].loopNext = newLabel;
                f1.name[newLabel].species = hidden;
                currLoop = newLabel;
                }
            }
        }
    return 0;
}


/**
 *Not currently used, but contains information on pseudo-potentials
*/
inta buildMetric( double latticePause,inta Z,   function_label func,inta am,   metric_label * metric){
    inta dim,dimB,nMet = 0,ai,space;
    const inta psu_stride = 17;
    const inta psu_length = 11;
    double lda[] = {
        1, 1, 0, 2, 2, 0.2, -4.06633, 0.677832, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        3, 3, 0, 4, 4, 0.4, -14.0094, 9.50991, -1.75327, 0.0834586, 0, 0, 0,0, 0, 0, 0,
        4, 4, 0, 4, 4, 0.325, -23.991, 17.1718, -3.31896,0.165083, 0, 0, 0, 0, 0, 0, 0,
        5, 3, 0, 3, 2, 0.4325, -5.60048,0.806284, 0, 0, 1, 0.373882, 6.23522, 0, 0, 0, 0,
        6, 4, 0, 3, 2,0.346473, -8.57533, 1.23413, 0, 0, 1, 0.304523, 9.53419, 0, 0, 0, 0,
        7, 5, 0, 3, 2, 0.288905, -12.2046, 1.75582, 0, 0, 1, 0.256912,0.256912, 0, 0, 0, 0,
        8, 6, 0, 3, 2, 0.247754, -16.4822, 2.37014, 0,0, 1, 0.222203, 18.1996, 0, 0, 0, 0,
        13, 3, 0, 6, 1, 0.45, -6.83406,0, 0, 0, 2, 0.465436, 2.81408, 1.93952, 3, 0.546243, 1.91601,
        14, 4,0, 6, 1, 0.44, -6.91363, 0, 0, 0, 2, 0.424334, 3.20813, 2.58888, 3,0.485359,2.65622,
        15, 5, 0, 6, 1, 0.43, -6.64097, 0, 0, 0, 2, 0.390738, 3.65826, 3.15066, 3, 0.440846, 3.28594,
        16, 6, 0, 6, 1,0.42, -6.59607, 0, 0, 0, 2, 0.362614, 4.22284, 3.66966, 3, 0.405311,3.88535};
    double blyp[] = {
        1, 1, 0, 2, 2, 0.2, -4.10561, 0.692787, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        3, 3, 0, 4, 4, 0.4, -14.1026, 9.65027, -1.79063, 0.0857313, 0, 0, 0,0, 0, 0, 0,
        4, 4, 0, 4, 4, 0.325, -24.0586, 17.2529, -3.33239,0.165305, 0, 0, 0, 0, 0, 0, 0,
        5, 3, 0, 3, 2, 0.424087, -6.08744,0.980916, 0, 0, 1, 0.371141, 6.32735, 0, 0, 0, 0,
        6, 4, 0, 3, 2,0.337633, -9.12847, 1.42513, 0, 0, 1, 0.302528, 9.65073, 0, 0, 0, 0,
        7, 5, 0, 3, 2, 0.281959, -12.7548, 1.94859, 0, 0, 1, 0.255444,13.6594, 0, 0, 0, 0,
        8, 6, 0, 3, 2, 0.24245, -17.0171, 2.56133, 0, 0,1, 0.221084, 18.3556, 0, 0, 0, 0,
        13, 3, 0, 6, 1, 0.45, -5.54822, 0,0, 0, 2, 0.505838, 3.02008, 1.06418, 3, 0.577572, 1.53528,
        14, 4, 0,6, 1, 0.44, -5.97966, 0, 0, 0, 2, 0.444927, 3.4402, 1.88129, 3,0.503637, 2.28821,
        15, 5, 0, 6, 1, 0.43, -5.87283, 0, 0, 0, 2,0.403545, 3.8762, 2.54131, 3, 0.452751, 2.9405,
        16, 6, 0, 6, 1, 0.42,-6.0083, 0, 0, 0, 2, 0.370401, 4.37362, 3.19573, 3, 0.413079, 3.5911};
    double *def, *PSU;
    
    
    switch ( func.fn ){
        case Pseudo:
            if ( latticePause > 1./func.param[2] ){
                if ( am < 1 )
                    exit(1);
                metric[nMet].fn.fn = Pseudo;
                metric[nMet].fn.param[0] = -func.param[0];
                metric[nMet].fn.param[2] = func.param[2];
                metric[nMet].fn.interval = func.interval;
                metric[nMet].fn.contr    = func.contr;
                metric[nMet].metric = interval;
                metric[nMet].beta[0] = 0;
                metric[nMet].beta[1] = 1./func.param[2];
                for ( space = 0; space < SPACE ; space++){
                    metric[nMet].pow[space] = 0;
                    metric[nMet].deriv[space] = 0;
                }
                nMet++;
            }else {
                if ( am < 2 )
                    exit(1);
                metric[nMet].fn.fn = Pseudo;
                metric[nMet].fn.param[0] = -func.param[0];
                metric[nMet].fn.param[2] = func.param[2];
                metric[nMet].fn.interval     = func.interval;
                metric[nMet].fn.contr    = func.contr;
                metric[nMet].metric = interval;
                metric[nMet].beta[0] = 0;
                metric[nMet].beta[1] = latticePause;
                for ( space = 0; space < SPACE ; space++){
                    metric[nMet].pow[space] = 0;
                    metric[nMet].deriv[space] = 0;
                }

                nMet++;
                metric[nMet].fn.fn = Pseudo;
                metric[nMet].fn.param[0] = -func.param[0];
                metric[nMet].fn.param[2] = func.param[2];
                metric[nMet].fn.interval = func.interval;
                metric[nMet].fn.contr    = func.contr;
                metric[nMet].metric = interval;
                metric[nMet].beta[0] = latticePause;
                metric[nMet].beta[1] = 1./func.param[2];
                for ( space = 0; space < SPACE ; space++){
                    metric[nMet].pow[space] = 0;
                    metric[nMet].deriv[space] = 0;
                }

                nMet++;
            }
            break;
        case Coulomb:
        case Yukawa:
            if ( am < 2 )
                exit(1);
            metric[nMet].fn.fn = func.fn;
            metric[nMet].fn.param[0] = -func.param[0];
            metric[nMet].fn.param[2] = func.param[2];
            metric[nMet].fn.interval = func.interval;
            metric[nMet].fn.contr    = func.contr;
            metric[nMet].metric = interval;
            metric[nMet].beta[0] = 0;
            metric[nMet].beta[1] = latticePause;
            for ( space = 0; space < SPACE ; space++){
                metric[nMet].pow[space] = 0;
                metric[nMet].deriv[space] = 0;
            }

            nMet++;
            metric[nMet].fn.fn = func.fn;
            metric[nMet].fn.param[0] = -func.param[0];
            metric[nMet].fn.param[2] = func.param[2];
            metric[nMet].fn.interval = func.interval;
            metric[nMet].fn.contr    = func.contr;
            metric[nMet].metric = semiIndefinite;
            metric[nMet].beta[0] = latticePause;
            for ( space = 0; space < SPACE ; space++){
                metric[nMet].pow[space] = 0;
                metric[nMet].deriv[space] = 0;
            }

            nMet++;
            break;
        case LDA:
        case BLYP:
           if ( LDA == func.fn )
                def = lda;
            else if  ( BLYP == func.fn )
                def = blyp;
            else {
                printf("not known \n");
                exit(0);
            }
            PSU = NULL;
            for ( ai = 0; ai < psu_length ; ai++){
                if ( def[psu_stride*ai] ==  abs(Z))//one place Z actually matters...pseudo Z.
                    PSU = def+psu_stride * ai;
            }
            if ( PSU == NULL )
            {
                printf ("not coded\n");
                exit(0);
            }
            //add pseudo core.
            if ( latticePause > 1./(sqrt(2.)*PSU[5]) ){
                if ( am < 1 )
                    exit(1);
                metric[nMet].fn.fn = Pseudo;
                metric[nMet].fn.param[0] = -PSU[1];
                metric[nMet].fn.param[2] = (sqrt(2.)*PSU[5]);
                metric[nMet].fn.interval = func.interval;
                metric[nMet].fn.contr    = func.contr;
                metric[nMet].metric = interval;
                metric[nMet].beta[0] = 0;
                metric[nMet].beta[1] =  1./(sqrt(2.)*PSU[5]);
                for ( space = 0; space < SPACE ; space++){
                    metric[nMet].pow[space] = 0;
                    metric[nMet].deriv[space] = 0;
                }

                nMet++;
            }else {
                if ( am < 2 )
                    exit(1);
                metric[nMet].fn.fn = Pseudo;
                metric[nMet].fn.param[0] = -PSU[1];
                metric[nMet].fn.param[2] = (sqrt(2.)*PSU[5]);
                metric[nMet].fn.interval = func.interval;
                metric[nMet].fn.contr    = func.contr;
                metric[nMet].metric = interval;
                metric[nMet].beta[0] = 0;
                metric[nMet].beta[1] = latticePause;
                for ( space = 0; space < SPACE ; space++){
                    metric[nMet].pow[space] = 0;
                    metric[nMet].deriv[space] = 0;
                }

                nMet++;
                metric[nMet].fn.fn = Pseudo;
                metric[nMet].fn.param[0] = -PSU[1];
                metric[nMet].fn.param[2] = (sqrt(2.)*PSU[5]);
                metric[nMet].fn.interval = func.interval;
                metric[nMet].fn.contr    = func.contr;
                metric[nMet].metric = interval;
                metric[nMet].beta[0] = latticePause;
                metric[nMet].beta[1] =  1./(sqrt(2.)*PSU[5]);
                for ( space = 0; space < SPACE ; space++){
                    metric[nMet].pow[space] = 0;
                    metric[nMet].deriv[space] = 0;
                }

                nMet++;
            }
            if (fabs(PSU[6])> 0. ){
                //single gaussian
                metric[nMet].fn.fn = Gaussian;
                metric[nMet].fn.param[0] = PSU[6]/pow(PSU[5], 0);
                metric[nMet].fn.interval = 0;
                metric[nMet].metric = dirac;
                metric[nMet].beta[0] =  1./(sqrt(2.)*PSU[5]);
                for ( space = 0; space < SPACE ; space++){
                    metric[nMet].pow[space] = 0;
                    metric[nMet].deriv[space] = 0;
                }

                nMet++;
                if ( nMet > am )
                    exit(1);
            }
            if (fabs(PSU[7])> 0. ){
                //3 gaussians.
                for ( dim = 0 ; dim < 3 ; dim++){
                    metric[nMet].fn.fn = Gaussian;
                    metric[nMet].fn.param[0] = PSU[7]/pow(PSU[5], 2);
                    metric[nMet].fn.interval = 0;
                    metric[nMet].metric = dirac;
                    metric[nMet].beta[0] =  1./(sqrt(2.)*PSU[5]);
                    for ( space = 0; space < SPACE ; space++){
                        metric[nMet].pow[space] = 0;
                        metric[nMet].deriv[space] = 0;
                    }
                    metric[nMet].pow[dim] = 2;
                    
                    nMet++;
                    if ( nMet > am )
                        exit(1);
                }
            }
            if (fabs(PSU[8])> 0. ){
                //6 gaussians
                for ( dim = 0 ; dim < 6 ; dim++){
                    metric[nMet].fn.fn = Gaussian;
                    metric[nMet].fn.param[0]= PSU[8]/pow(PSU[5], 4);
                    metric[nMet].fn.interval = 0;
                    metric[nMet].metric = dirac;
                    metric[nMet].beta[0] =  1./(sqrt(2.)*PSU[5]);
                    for ( space = 0; space < SPACE ; space++){
                        metric[nMet].pow[space] = 0;
                        metric[nMet].deriv[space] = 0;
                    }
                    //a^2 + 2 a b + b^2 + 2 a c + 2 b c + c^2

                    switch( dim ){
                        case 0:
                            metric[nMet].pow[0] = 4;
                            break;
                        case 1:
                            metric[nMet].fn.param[0] *= 2;
                            metric[nMet].pow[0] = 2;
                            metric[nMet].pow[1] = 2;
                            break;
                        case 2:
                            metric[nMet].pow[1] = 4;
                            break;
                        case 3:
                            metric[nMet].fn.param[0] *= 2;
                            metric[nMet].pow[1] = 2;
                            metric[nMet].pow[2] = 2;
                            break;
                        case 4:
                            metric[nMet].pow[2] =4;
                            break;
                        case 5:
                            metric[nMet].fn.param[0] *= 2;
                            metric[nMet].pow[0] = 2;
                            metric[nMet].pow[2] = 2;
                            break;
                            
                    }
                    
                    nMet++;
                    if ( nMet > am )
                        exit(1);

                }
            }
            if (fabs(PSU[9])> 0. ){
                //6 gaussians
                if ( SPACE < 3 ){
                    printf("ack!");
                    exit(1);
                }
                for ( dim = 0 ; dim < 10 ; dim++){
                    metric[nMet].fn.fn = Gaussian;
                    metric[nMet].fn.param[0] = PSU[9]/pow(PSU[5], 6);
                    metric[nMet].fn.interval = 0;
                    metric[nMet].metric = dirac;
                    metric[nMet].beta[0] =  1./(sqrt(2.)*PSU[5]);
                    for ( space = 0; space < SPACE ; space++){
                        metric[nMet].pow[space] = 0;
                        metric[nMet].deriv[space] = 0;
                    }
                    //a^3 + 3 a^2 b + 3 a b^2 + b^3 + 3 a^2 c + 6 a b c + 3 b^2 c +
                    //3 a c^2 + 3 b c^2 + c^3
                    

                    switch( dim ){
                        case 0:
                            metric[nMet].pow[0] = 6;
                            break;
                        case 1:
                            metric[nMet].fn.param[0] *= 3;
                            metric[nMet].pow[0] = 4;
                            metric[nMet].pow[1] = 2;
                            break;
                        case 2:
                            metric[nMet].fn.param[0] *= 3;
                            metric[nMet].pow[0] = 2;
                            metric[nMet].pow[1] = 4;
                            break;
                        case 3:
                            metric[nMet].pow[1] = 6;
                            break;
                        case 4:
                            metric[nMet].fn.param[0] *= 3;
                            metric[nMet].pow[0] = 4;
                            metric[nMet].pow[2] = 2;
                            break;
                        case 5:
                            metric[nMet].fn.param[0] *= 6;
                            metric[nMet].pow[0] = 2;
                            metric[nMet].pow[1] = 2;
                            metric[nMet].pow[2] = 2;
                            break;
                        case 6:
                            metric[nMet].fn.param[0] *= 3;
                            metric[nMet].pow[1] = 4;
                            metric[nMet].pow[2] = 2;
                            break;
                        case 7:
                            metric[nMet].fn.param[0] *= 3;
                            metric[nMet].pow[0] = 2;
                            metric[nMet].pow[2] = 4;
                            break;
                        case 8:
                            metric[nMet].fn.param[0] *= 3;
                            metric[nMet].pow[1] = 2;
                            metric[nMet].pow[2] = 4;
                            break;
                        case 9:
                            metric[nMet].pow[2] = 6;
                            break;
                    }
                
                    nMet++;
                    if ( nMet > am )
                        exit(1);
                    
                }
            }
            if (fabs(PSU[10]) > 0. ){
                //single separated gaussian
                metric[nMet].fn.fn = Gaussian;
                metric[nMet].fn.param[0] = PSU[12]/pow(sqrt(pi)*PSU[11],3);
                metric[nMet].fn.interval= 0;
                metric[nMet].metric = separateDirac;
                metric[nMet].beta[0] =  1./(sqrt(2.)*PSU[5]);
                for ( space = 0; space < SPACE ; space++){
                    metric[nMet].pow[space] = 0;
                    metric[nMet].powB[space] = 0;
                    metric[nMet].deriv[space] = 0;
                }

                nMet++;
                if ( nMet > am )
                    exit(1);
            }
            if (fabs(PSU[10]) > 1. ){
                //9 separated radial gaussian
                for ( dim = 0 ;dim < 3 ;dim++)
                    for ( dimB = 0 ;dimB < 3 ;dimB++){
                        metric[nMet].fn.fn = Gaussian;
                        metric[nMet].fn.param[0] = 2./15.*PSU[16]/pow(pi,1.500)/pow(PSU[15],7.);
                        metric[nMet].fn.interval = 0;
                        metric[nMet].metric = separateDirac;
                        metric[nMet].beta[0] =  1./(sqrt(2.)*PSU[11]);
                        for ( space = 0; space < SPACE ; space++){
                            metric[nMet].pow[space] = 0;
                            metric[nMet].powB[space] = 0;
                            metric[nMet].deriv[space] = 0;
                        }
                        //(x2 + y2 + z2 ) ( x'2 + y'2 + z'2 )
                        metric[nMet].pow[dim] = 2;
                        metric[nMet].powB[dimB] = 2;
                        
                        nMet++;
                        if ( nMet > am )
                            exit(1);
                    }
            }
            if (fabs(PSU[14]) > 0. ){
                //3 separated p-gaussian
                for ( dim = 0 ;dim < 3 ;dim++){
                    metric[nMet].fn.fn = Gaussian;
                    metric[nMet].fn.param[0] = 2.*PSU[16]/pow(pi,1.500)/pow(PSU[15],5.);
                    metric[nMet].fn.interval = 0;
                    metric[nMet].metric = separateDirac;
                    metric[nMet].beta[0] =  1./(sqrt(2.)*PSU[15]);
                    for ( space = 0; space < SPACE ; space++){
                        metric[nMet].pow[space] = 0;
                        metric[nMet].deriv[space] = 0;
                    }
                    metric[nMet].pow[dim]  = 1;
                    metric[nMet].powB[dim] = 1;

                    nMet++;
                    if ( nMet > am )
                        exit(1);
                }
            }
            break;
        case Gaussian:
            
                //single gaussian
                metric[nMet].fn.fn = Gaussian;
                metric[nMet].fn.param[0] = func.param[0];
            metric[nMet].fn.param[2] = func.param[2];
                metric[nMet].fn.interval = 0;
                metric[nMet].metric = dirac;
                metric[nMet].beta[0] = func.param[2];
                for ( space = 0; space < SPACE ; space++){
                    metric[nMet].pow[space] = 0;
                    metric[nMet].deriv[space] = 0;
                }
                nMet++;
                if ( nMet > am )
                    exit(1);
            
            break;
        default:
            break;

    }
    
    return nMet;
}


/**
 *Eikon of Kinetic
 *
 * MAY NEED A BLOCH K
 *
 *@param c1 parameters
 *@param f1 container
 *@param scalar overall scalar multiply
 *@param invert switch to turn on a particle-1 inversion
 *@param act Symmetry Adaption related work, for group action
 *@param bl the address of the interaction, i.e. particle-1 or particle-12
 *@param[in] single linked list
 *@param label the component group ID
 *@param overline legacy
 *@param cmpl make it real for now
 */
inta buildKinetic( calculation *c1, sinc_label *f1,double scalar,inta invert,inta act, blockType bl, division single,inta label, inta overline, spinType cmpl){
    inta dim,spacy,id = 0;
      division li = single;
    while ( f1->name[li].chainNext != nullName)
        li =f1->name[li].chainNext;
    //x is beta.
    if ( overline ){
        printf("periodic boundary conditions not implemented yet!");
        exit(0);
    }
    printf("Kinetic\tx %d act %d block %d (%f)\n", label,act,bl,scalar);

      division headLabel,currLabel;
    currLabel = li;
    for ( dim = 0 ; dim < SPACE ; dim++)
        if ( f1->canon[dim].body != nada)
        {
            if ( f1->canon[dim].label == label)
                {
                    
                    headLabel = anotherLabel(f1,0,nada);
                    f1->name[currLabel].chainNext = headLabel;
                    f1->name[headLabel].species = eikon;
                    //new term
                      division memoryLabel = anotherLabel(f1,all,one);
                    f1->name[headLabel].loopNext = memoryLabel;
                    f1->name[memoryLabel].species = eikonKinetic;
                    f1->name[memoryLabel].Current[0] = 1;
                    f1->name[headLabel].multId = id++;//aloways unique

            
                    for (spacy = 0 ; spacy < SPACE ; spacy++)//set term across basis
                        if ( f1->canon[spacy].body != nada){
                            f1->name[memoryLabel].space[spacy].act = act;
                            if ( f1->canon[spacy].label == label && spacy == dim)
                                {
                                    f1->name[memoryLabel].space[spacy].body = one;
                                    f1->name[memoryLabel].space[spacy].block = bl;
                                    streams(*f1, memoryLabel, 0, spacy)[0] = -0.500*scalar;
                                }else{
                                        f1->name[memoryLabel].space[spacy].block = id0;
                                }
                    
                        }
                    currLabel = headLabel;
            }
        }
        
        
        
        
    return 0;
}

/**
 *Eikon of a constant function
 *
 *@param c1 parameters
 *@param f1 container
 *@param scalar overall scalar multiply
 *@param invert switch to turn on a particle-1 inversion
 *@param act Symmetry Adaption related work, for group action
 *@param bl the address of the interaction, i.e. particle-1 or particle-12
 *@param[in] single linked list
 *@param label the component group ID
 *@param overline legacy
 *@param cmpl make it real for now
*/
inta buildConstant(  calculation *c1,   sinc_label *f1,double scalar,inta invert,inta act,   blockType bl,  division single,inta label, inta overline,   spinType cmpl){
    inta space,id=0;
      division li = single;
    while ( f1->name[li].chainNext != nullName)
        li =f1->name[li].chainNext;
    //x is beta.
    if ( overline ){
        printf("periodic boundary conditions not implemented yet!");
        exit(0);
    }
    printf("constant %f invert %d act %d block %d\n", scalar,invert,act,bl);
      division headLabel;
    headLabel = anotherLabel(f1,0,nada);
    f1->name[li].chainNext = headLabel;
    f1->name[headLabel].species = eikon;
      division memoryLabel = anotherLabel(f1,all,one);
    f1->name[headLabel].loopNext = memoryLabel;
    f1->name[memoryLabel].species = eikonConstant;
    f1->name[memoryLabel].Current[0] = 1;
    f1->name[memoryLabel].Partition = 1;
    f1->name[headLabel].multId = id++;

    for ( space = 0 ; space < SPACE ; space++)
        if ( f1->canon[space].body != nada){
            if ( f1->canon[space].label == label)
                {
                    f1->name[memoryLabel].space[space].block = bl;
                    f1->name[memoryLabel].space[space].act = act;
                    streams(*f1, memoryLabel, 0, space)[0] = scalar;//redundant and irrelevant.
                }else
                {
                    f1->name[memoryLabel].space[space].block = id0;
                    f1->name[memoryLabel].space[space].act = 1;
                }
            }
    return 0;
}


/**
 *Eikon of linear X function
 *
 *@param c1 parameters
 *@param f1 container
 *@param scalar overall scalar multiply
 *@param invert switch to turn on a particle-1 inversion
 *@param act Symmetry Adaption related work, for group action
 *@param bl the address of the interaction, i.e. particle-1 or particle-12
 *@param[in] single linked list
 *@param label the component group ID
 *@param overline legacy
 *@param cmpl make it real for now
*/
inta buildLinear(  calculation *c1,   sinc_label *f1,double scalar,inta invert,inta act,   blockType bl,  division single,inta label, inta overline,   spinType cmpl){
    inta dim,spacy,id=0;
      division li = single;
    while ( f1->name[li].chainNext != nullName)
        li =f1->name[li].chainNext;
    //x is beta.
    if ( overline ){
        printf("periodic boundary conditions not implemented yet!");
        exit(0);
    }
    printf("Force\tx %d act %d block %d \t(%f)\n", label,act,bl,scalar);

         division headLabel,currLabel;
       currLabel = li;
       for ( dim = 0 ; dim < SPACE ; dim++)
           if ( f1->canon[dim].body != nada)
           {
               if ( f1->canon[dim].label == label)
                   {
                       
                       headLabel = anotherLabel(f1,0,nada);
                       f1->name[currLabel].chainNext = headLabel;
                       f1->name[headLabel].species = eikon;
                       f1->name[headLabel].multId = id++;
                       //new term
                         division memoryLabel = anotherLabel(f1,all,one);
                       f1->name[headLabel].loopNext = memoryLabel;
                       f1->name[memoryLabel].species = eikonLinear;
                       f1->name[memoryLabel].Current[0] = 1;

               
                       for (spacy = 0 ; spacy < SPACE ; spacy++)//set term across basis
                           if ( f1->canon[spacy].body != nada){
                               f1->name[memoryLabel].space[spacy].act = 1;//i think this is weak, but not concerning now...
                               if ( f1->canon[spacy].label == label && spacy == dim)
                                   {
                                       f1->name[memoryLabel].space[spacy].body = one;
                                       f1->name[memoryLabel].space[spacy].block = bl;
                                       streams(*f1, memoryLabel, 0, spacy)[0] = scalar;
                                   }else{
                                           f1->name[memoryLabel].space[spacy].block = id0;
                                   }
                       
                           }
                       currLabel = headLabel;
               }
           }

    return 0;
}


/**
 *Eikon of deriv operator
 *
 *@param c1 parameters
 *@param f1 container
 *@param scalar overall scalar multiply
 *@param invert switch to turn on a particle-1 inversion
 *@param act Symmetry Adaption related work, for group action
 *@param bl the address of the interaction, i.e. particle-1 or particle-12
 *@param[in] single linked list
 *@param label the component group ID
 *@param overline legacy
 *@param cmpl make it real for now
*/
inta buildDeriv(  calculation *c1,   sinc_label *f1,double scalar,inta invert,inta act,   blockType bl,  division single,inta label, inta overline,   spinType cmpl){
    inta dim,spacy,id=0;
      division li = single;
    while ( f1->name[li].chainNext != nullName)
        li =f1->name[li].chainNext;
    //x is beta.
    if ( overline ){
        printf("periodic boundary conditions not implemented yet!");
        exit(0);
    }
    printf("Force\tx %d act %d block %d \t(%f)\n", label,act,bl,scalar);

         division headLabel,currLabel;
       currLabel = li;
       for ( dim = 0 ; dim < SPACE ; dim++)
           if ( f1->canon[dim].body != nada)
           {
               if ( f1->canon[dim].label == label)
                   {
                       
                       headLabel = anotherLabel(f1,0,nada);
                       f1->name[currLabel].chainNext = headLabel;
                       f1->name[headLabel].species = eikon;
                       f1->name[headLabel].multId = id++;
                       //new term
                         division memoryLabel = anotherLabel(f1,all,one);
                       f1->name[headLabel].loopNext = memoryLabel;
                       f1->name[memoryLabel].species = eikonDeriv;
                       f1->name[memoryLabel].Current[0] = 1;

                       ///set term across basis
                       for (spacy = 0 ; spacy < SPACE ; spacy++)
                           if ( f1->canon[spacy].body != nada){
                               f1->name[memoryLabel].space[spacy].act = 1;//i think this is weak, but not concerning now...
                               if ( f1->canon[spacy].label == label && spacy == dim)
                                   {
                                       f1->name[memoryLabel].space[spacy].body = one;
                                       f1->name[memoryLabel].space[spacy].block = bl;
                                       streams(*f1, memoryLabel, 0, spacy)[0] = scalar;
                                   }else{
                                           f1->name[memoryLabel].space[spacy].block = id0;
                                   }
                       
                           }
                       currLabel = headLabel;
               }
           }

    return 0;
}


/**
 *Eikon of state elements
 *
 *@param c1 parameters
 *@param f1 container
 *@param scalar overall scalar multiply
 *@param invert switch to turn on a particle-1 inversion
 *@param act Symmetry Adaption related work, for group action
 *@param bl the address of the interaction, i.e. particle-1 or particle-12
 *@param[in] single linked list
 *@param label the component group ID
 *@param overline legacy
 *@param cmpl make it real for now
 *@param bra an index of an element
 *@param ket an index of an element
*/
inta buildElement(  calculation *c1,   sinc_label *f1,double scalar,inta invert,inta act,   blockType bl,  division single,inta label, inta overline,   spinType cmpl,inta bra, inta ket){
    inta space;
      division li = single;
    while ( f1->name[li].chainNext != nullName)
        li =f1->name[li].chainNext;
    //x is beta.
    if ( overline ){
        printf("periodic boundary conditions not implemented yet!");
        exit(0);
    }
    printf("Element\tx %d act %d block %d \t(%d,%d)->(%f)\n", label,act,bl,bra,ket,scalar);

      division headLabel;
    headLabel = anotherLabel(f1,0,nada);
    f1->name[li].chainNext = headLabel;
    f1->name[headLabel].species = eikon;
      division memoryLabel = anotherLabel(f1,all,one);
    f1->name[headLabel].loopNext = memoryLabel;
    f1->name[headLabel].multId = 0;

    f1->name[memoryLabel].species = eikonElement;
    f1->name[memoryLabel].Current[0] = 1;

    for ( space = 0 ; space < SPACE ; space++)
        if ( f1->canon[space].body != nada){
            f1->name[headLabel].space[space].act = act;
            if ( f1->canon[space].label == label)
                {
                    f1->name[memoryLabel].space[space].body = one;

                    streams(*f1, memoryLabel, 0, space)[0] = scalar;
                    f1->name[memoryLabel].space[space].block = tv1;
                    f1->name[memoryLabel].space[space].bra = bra;
                    f1->name[memoryLabel].space[space].ket = ket;
                }else
                {
                    f1->name[memoryLabel].space[space].block = id0;
                }
            }

    return 0;
}




/**
 *Eikon of square X function
 *
 *@param c1 parameters
 *@param f1 container
 *@param scalar overall scalar multiply
 *@param invert switch to turn on a particle-1 inversion
 *@param act Symmetry Adaption related work, for group action
 *@param bl the address of the interaction, i.e. particle-1 or particle-12
 *@param[in] single linked list
 *@param label the component group ID
 *@param overline legacy
 *@param cmpl make it real for now
*/
inta buildSpring(  calculation *c1,   sinc_label *f1,double scalar,inta invert,inta act,   blockType bl,  division single,inta label, inta overline,   spinType cmpl){
    inta dim,spacy,id=0;
      division li = single;
    while ( f1->name[li].chainNext != nullName)
        li =f1->name[li].chainNext;
    //x is beta.
    if ( overline ){
        printf("periodic boundary conditions not implemented yet!");
        exit(0);
    }
    printf("Spring\tx %d act %d block %d (%f)\n", label,act,bl,scalar);

    
         division headLabel,currLabel;
       currLabel = li;
       for ( dim = 0 ; dim < SPACE ; dim++)
           if ( f1->canon[dim].body != nada)
           {
               if ( f1->canon[dim].label == label)
                   {

                       headLabel = anotherLabel(f1,0,nada);
                       f1->name[currLabel].chainNext = headLabel;
                       f1->name[headLabel].species = eikon;
                       //new term
                       f1->name[headLabel].multId = id++;

                         division memoryLabel = anotherLabel(f1,all,one);
                       f1->name[headLabel].loopNext = memoryLabel;
                       f1->name[memoryLabel].species = eikonSpring;
                       f1->name[memoryLabel].Current[0] = 1;

               
                       for (spacy = 0 ; spacy < SPACE ; spacy++)//set term across basis
                           if ( f1->canon[spacy].body != nada){
                               f1->name[memoryLabel].space[spacy].act = act;
                               if ( f1->canon[spacy].label == label && spacy == dim)
                                   {
                                       f1->name[memoryLabel].space[spacy].body = one;
                                       f1->name[memoryLabel].space[spacy].block = bl;
                                       streams(*f1, memoryLabel, 0, spacy)[0] = 0.500*scalar;
                                   }else{
                                           f1->name[memoryLabel].space[spacy].block = id0;
                                   }
                       
                           }
                       currLabel = headLabel;
               }
           }
               return 0;
}


/**
 *Eikon of oneBody interaction
 *
 *@param c1 parameters
 *@param f1 container
 *@param scalar overall scalar multiply
 *@param invert switch to turn on a particle-1 inversion
 *@param act Symmetry Adaption related work, for group action
 *@param bl the address of the interaction, i.e. particle-1 or particle-12
 *@param adjustOne to change the defined basis lattice
 *@param[in] single linked list
 *@param particle1 the component group ID
 *@param embed the number of trivial argumented Gaussians
 *@param overline legacy
 *@param cmpl make it real for now
 *@param mu the metric
 *@param a geometry/Z of a-th ion
 */
inta buildExternalPotential(  calculation *c1,   sinc_label *f1,double scalar, inta invert,inta act,  blockType bl, double adjustOne,   division single,  inta particle1,inta embed, inta overline,   spinType cmpl,  metric_label mu,inta a){
    inta ra=0;
    printf("oneBody act %d block %d atom %d >%d<-- (%f)\n", act,bl, a,embed,scalar);

            if ( mu.metric == interval || mu.metric == semiIndefinite)
                ra += (mu.fn.interval);
            else if ( mu.metric == dirac )
                ra++;
        if ( bootedQ(*f1) ){
                    separateInteraction(f1,scalar*c1->i.atoms[a].Z, c1->i.atoms[a].position+1,invert,act,bl, adjustOne, single, mu, cmpl,  0, particle1,one,embed);
                }
        
    if ( bootedQ(*f1) ){
    }

    return ra;
}

/**
 *Eikon of twoBody interaction
 *
 *@param c1 parameters
 *@param f1 container
 *@param scalar overall scalar multiply
 *@param invert switch to turn on a particle-1 inversion
 *@param act Symmetry Adaption related work, for group action
 *@param bl the address of the interaction, i.e. particle-1 or particle-12
 *@param adjustOne to change the defined basis lattice
 *@param[in] pair linked list
 *@param particle1 the component group ID
 *@param embed the number of trivial argumented Gaussians
 *@param overline legacy
 *@param cmpl make it real for now
 *@param mu the metric
*/
inta buildPairWisePotential(  calculation *c1,   sinc_label *f1,double scalar,inta invert, inta act,  blockType bl, double adjustOne,  division pair,  inta particle1 ,inta embed, inta overline,   spinType cmpl,  metric_label mu){
    inta ra=0;
    printf("twoBody act %d block %d [adjust (%f)] >%d<-- (%f)\n",act,bl,adjustOne,embed,scalar);
    
    if ( mu.metric == interval || mu.metric == semiIndefinite)
        ra += (mu.fn.interval);
    else if ( mu.metric == dirac )
        ra++;

    if ( bootedQ(*f1)  ){
            double zero[6];
            zero[0] = 0.;
            zero[1] = 0.;
            zero[2] = 0.;
            zero[3] = 0.;
            zero[4] = 0.;
            zero[5] = 0.;
            separateInteraction(f1, scalar,zero,invert,act,bl, adjustOne, pair , mu, cmpl, 0, particle1,two,embed);
        }
    return ra;
}


/**
 *Diagonal Matrix reference
 *
 *@param f field
 *@param filename file to put into vectorDiagonalMatrix, be sure to switch it on!
 *@param[in] single linked list
 */
inta assignDiagonalMatrix(calculation *c1,   field *f, char * filename, division single){
    sinc_label * f1 = &f->f;
    inta index = 0;
    division li = single;
    while ( f1->name[li].chainNext != nullName)
        li =f1->name[li].chainNext;
    printf("diagonalMatrix %s\n",filename);
    division headLabel;
    headLabel = anotherLabel(f1,0,nada);
    f1->name[li].chainNext = headLabel;
    f1->name[headLabel].Current[0] = 1;
 
    f1->name[headLabel].loopNext = vectorDiagonalMatrix;
    f1->name[headLabel].multId = 101;
    
    f1->name[vectorDiagonalMatrix].species = diagonalMatrix;
#ifndef APPLE
    tLoadEigenWeights(c1,*f,filename, &index, vectorDiagonalMatrix, 0);
#else
    tBoot(f->f, vectorDiagonalMatrix, 0, 1);
#endif
    return 0;
}
