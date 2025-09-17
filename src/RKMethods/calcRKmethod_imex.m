function method = calcRKmethod_imex(MethodName)

switch upper(MethodName)
    
    case 'EULER'
        
        method.aI   = [ 0.0     0.0 
                        0.0     1.0 ];
        method.c    = [ 0.0;    1.0 ];
        method.s            = 2;
        method.GSA          = true;
        method.const_diag   = true;
        method.order        = 1;
        
    case 'ARS443'
    
        method.aI   = [ 0.0 	0.0 	0.0 	0.0 	0.0 ;
						0.0 	1/2 	0.0 	0.0 	0.0 ;
						0.0 	1/6 	1/2 	0.0 	0.0 ;
						0.0 	-1/2 	1/2 	1/2 	0.0 ;
						0.0 	3/2 	-3/2 	1/2 	1/2 ];
		method.aE   = [ 0.0 	0.0 	0.0 	0.0  	0.0;
						1/2 	0.0 	0.0 	0.0  	0.0;
						11/18 	1/18 	0.0 	0.0  	0.0;
						5/6 	-5/6 	1/2 	0.0  	0.0;
                        1/4 	7/4 	3/4 	-7/4  	0.0 ];
        method.c        = sum(method.aI, 2);
        method.s        = length(method.c);
        method.bI       = method.aI(end,:).';
        method.bE       = method.aE(end,:).';
        method.bhatI    = [ 0.0     1.0     0.0     0.0     0.0 ].';
        method.bhatE    = [ 0.0     1.0     0.0     0.0     0.0 ].';
        method.GSA      = true;
        method.const_diag   = true;
        method.order    = 3;
        
    case 'BPR3'
        
        %3rd order RK BPR(3,5,3) from Boscarino et al (2012):
        method.aI       = [ 0.0     0.0     0.0     0.0     0.0
                            0.5     0.5     0.0     0.0     0.0 
                            5/18    -1/9    0.5     0.0     0.0 
                            0.5     0.0     0.0     0.5     0.0 
                            0.25    0.0     0.75    -0.5    0.5 ];
        method.aE       = [ 0.0     0.0     0.0     0.0     0.0
                            1.0     0.0     0.0     0.0     0.0 
                            4/9     2/9     0.0     0.0     0.0 
                            0.25    0.0     0.75    0.0     0.0 
                            0.25    0.0     0.75    0.0     0.0 ];
        method.bI       = method.aI(end,:).';
        method.bE       = method.aE(end,:).';
        
        if false && t_evap<1
        % SHORT SIMULATION
            method.bhatI    = [ 0.0    1.0    0.0     0.0     0.0 ].';
            method.bhatE    = [ 0.0     1.0     0.0     0.0     0.0 ].';
        else
        % LONG SIMULATION (+-100s)
            method.bhatI    = [ 5/18    -1/9    0.5     0.0     0.0 ].';
            method.bhatE    = [ 4/9     2/9     0.0     0.0     0.0 ].';
        end
  
        method.GSA          = true;
        method.const_diag   = true;
        method.c            = sum(method.aI,2);
        method.s            = length(method.c);
        method.order        = 3;
        
    case 'KC4B'
        
        %Implicit method:
        gamma       = 1235/10000;
        aI          = zeros(8,8);
        aI(3,2)     = 624185399699/4186980696204;
        aI(4,2)     = 1258591069120/10082082980243;
        aI(4,3)     = -322722984531/8455138723562;
        aI(5,2)     = -436103496990/5971407786587;
        aI(5,3)     = -2689175662187/11046760208243;
        aI(5,4)     = 4431412449334/12995360898505;
        aI(6,2)     = -2207373168298/14430576638973;
        aI(6,3)     = 242511121179/3358618340039;
        aI(6,4)     = 3145666661981/7780404714551;
        aI(6,5)     = 5882073923981/14490790706663;
        aI(7,2)     = 0.0;
        aI(7,3)     = 9164257142617/17756377923965;
        aI(7,4)     = -10812980402763/74029279521829;
        aI(7,5)     = 1335994250573/5691609445217;
        aI(7,6)     = 2273837961795/8368240463276;
        for kk=2:7
            aI(kk,kk)   = gamma;
            aI(kk,1)    = aI(kk,2);
        end
        aI(8,:)     = aI(7,:);

        %Explicit method:
        aE          = zeros(8,8);
        aE(2,1)     = 247/1000;
        aE(3,1)     = 247/4000;
        aE(3,2)     = 2694949928731/7487940209513;
        aE(4,1)     = 464650059369/8764239774964;
        aE(4,2)     = 878889893998/2444806327765;
        aE(4,3)     = -952945855348/12294611323341;
        aE(5,1)     = 476636172619/8159180917465;
        aE(5,2)     = -1271469283451/7793814740893;
        aE(5,3)     = -859560642026/4356155882851;
        aE(5,4)     = 1723805262919/4571918432560;
        aE(6,1)     = 6338158500785/11769362343261;
        aE(6,2)     = -4970555480458/10924838743837;
        aE(6,3)     = 3326578051521/2647936831840;
        aE(6,4)     = -880713585975/1841400956686;
        aE(6,5)     = -1428733748635/8843423958496;
        aE(7,1)     = 760814592956/3276306540349;
        aE(7,2)     = 760814592956/3276306540349;
        aE(7,3)     = -47223648122716/6934462133451;
        aE(7,4)     = 71187472546993/9669769126921;
        aE(7,5)     = -13330509492149/9695768672337;
        aE(7,6)     = 11565764226357/8513123442827;
        aE(8,:)     = aI(8,:);
    
        method.aI   = aI(1:end-1, 1:end-1);
        method.aE   = aE(1:end-1, 1:end-1);
        method.c    = sum(method.aI, 2);
        method.s    = length(method.c);
        method.bI   = aI(end,1:end-1).';
        method.bE   = aE(end,1:end-1).';
        method.GSA  = false;
        
    case 'KC35'
        
        Am          = zeros(5,5);
        sq2         = sqrt(2.0);
        for ii=2:5
            Am(ii,ii)   = 9/40;
        end
        Am(2,1)     = 9/40;
        Am(3,1)     = 9*(1+sq2)/80;
        Am(3,2)     = 9*(1+sq2)/80;
        Am(4,1)     = (22+15*sq2)/(80*(1+sq2));
        Am(4,2)     = (22+15*sq2)/(80*(1+sq2));
        Am(4,3)     = -7/(40*(1+sq2));
        Am(5,1)     = (2398+1205*sq2)/(2835*(4+3*sq2));
        Am(5,2)     = (2398+1205*sq2)/(2835*(4+3*sq2));
        Am(5,3)     = -2374*(1+2*sq2)/(2835*(5+3*sq2));
        Am(5,4)     = 5827/7560;

        bhat        = zeros(5,1);
        bhat(1)     = 4555948517383/24713416420891;
        bhat(2)     = 4555948517383/24713416420891;
        bhat(3)     = -7107561914881/25547637784726;
        bhat(4)     = 30698249/44052120;
        bhat(5)     = 49563/233080;
        
        method.GSA          = true;
        method.const_diag   = true;
        method.aI           = Am;
        method.bI           = Am(end,:).';
        method.bhatI        = bhat;
        method.c            = sum(method.aI,2);
        method.s            = length(method.c);
        method.order        = 3;

    otherwise
        error(['Unknown method ', MethodName])
end
end

