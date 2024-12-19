
function africa_pshape=coastline_africa(flag)

%--------------------------------------------------------------------------
% Polyshape object that is a coarse approximation of Africa (including
% Madagascar).
%--------------------------------------------------------------------------

% In vertices we store, respectively, Longitude and Latitute, in degrees.
% Column 1: Longitude. Column 2: Latitude.

madagascar_coord=[
    44.692595702572078 -16.172241632787330
    44.951115953285054 -16.124276543967774
    45.254998848595754 -15.946247665520591
    45.323783032117966 -15.966051316717168
    45.460792198416001 -16.025168564995823
    45.550881389007252 -15.994872826977952
    45.658566529235067 -15.835824862803094
    46.164865201239429 -15.723812976561669
    46.256132443170060 -15.944434052220792
    46.403001775808015 -15.941520309375202
    46.448520773489122 -15.849218319834041
    46.364568614866499 -15.630707105866994
    46.724329704667966 -15.401189940291488
    46.922366632418033 -15.329542146313946
    47.067652218311494 -15.432380173567227
    47.246264396131032 -15.329572247194939
    47.109774739565715 -15.100100014870064
    47.217996473668485 -14.947386177170534
    47.472718464919183 -14.708918582357889
    47.444391385160209 -14.906570489101453
    47.545743454582485 -14.931089312460454
    47.709287327118666 -14.655790509424977
    47.698823056768035 -14.580561034916274
    47.693907975783205 -14.425685893218573
    47.738942629324228 -14.293584661389083
    47.829274229125396 -14.296923795294257
    47.913286261455710 -14.257001401069839
    48.013094921669996 -14.111990702927583
    47.920997430471317 -13.725411384250929
    47.951551966606253 -13.574288687389124
    48.054054730633794 -13.552213557937328
    48.085662749977608 -13.558314729460553
    48.190300058988754 -13.778395314194000
    48.300625186431674 -13.819776617388014
    48.334924357826623 -13.537517226359050
    48.691560830818467 -13.405824574061686
    48.816995479760642 -13.297647543709006
    48.883165779054757 -13.114885597644031
    48.965172246868669 -12.839238790931347
    48.874405733118429 -12.642874401747568
    48.862793223576041 -12.520623864700369
    49.192215795886526 -12.198411715661068
    49.447235628484378 -12.362237340144906
    49.607753101584358 -12.675559862890795
    49.837774048679286 -12.836759530642079
    49.877551706897762 -12.934474465431780
    50.021182691999677 -13.318697418242648
    50.167802081358140 -14.115328168297115
    50.186796538938658 -14.376718159990951
    50.251447008131819 -14.751379476028154
    50.426662589488657 -15.087295980637192
    50.499448936094318 -15.302036905071610
    50.481409562100119 -15.473857932569995
    50.449556898891096 -15.547385766339556
    50.423324160388454 -15.594382764771911
    50.176701492296878 -15.977291460851598
    50.089412900207030 -15.913403783064654
    50.022580624103732 -15.757113965934504
    49.813865562251706 -15.424474084623581
    49.661119813958742 -15.472468265353635
    49.626025234022052 -15.637923169896014
    49.704816134370489 -15.994664514738872
    49.839879310091753 -16.265448428403651
    49.851663771441686 -16.455160582628977
    49.753897355921076 -16.716329050130401
    49.495453527221450 -17.121553943621496
    49.484175538736622 -17.320630476711127
    49.519789607878899 -17.711807510240149
    49.433480832276203 -18.142244697964866
    49.103939783225407 -19.105507363082925
    48.719367193310866 -20.274333071958782
    48.137749839306203 -21.998570992840573
    47.509935981337811 -24.029548775802066
    47.242984144637148 -24.735609859654680
    47.110037401936530 -24.852023404645745
    46.765213757568418 -25.142912853183365
    46.672950185639415 -25.192091096226722
    46.235357620507664 -25.243258516608165
    45.837793058501568 -25.393920693304558
    45.603915427913229 -25.548419912589562
    45.159174893173912 -25.596817259897378
    44.898585464776062 -25.392808498337665
    44.802500031788036 -25.332955338029706
    44.604311769803857 -25.259528330652401
    44.492528720999566 -25.295676555573152
    44.395452327797265 -25.274273149242479
    44.001190690910484 -24.963182549642944
    43.950137639017008 -24.859873853973241
    43.932209593335052 -24.622375702197420
    43.893812026719417 -24.540654413308964
    43.730972730929160 -24.452379921043217
    43.697020338457030 -24.372198988157404
    43.657344356341490 -24.044885210455380
    43.748709817153141 -23.548200292212968
    43.707048789787251 -23.387064040800997
    43.579893566640138 -23.145566417980120
    43.305513166907701 -22.675756886351113
    43.222090773884986 -22.428899381627875
    43.256100214683890 -22.102808301870919
    43.431892805733959 -21.433688268257296
    43.533108169934700 -21.332920242124498
    43.749762140456411 -21.255609047089916
    44.095914131827804 -20.658996814647079
    44.405683327748271 -20.148781453995952
    44.452629342908523 -20.006532873905176
    44.480101668255685 -19.521500738510490
    44.421380825808527 -19.358808367173364
    44.346449627513415 -19.257086172355169
    44.242205065345559 -19.090671188527327
    44.240533750198395 -19.005717431077528
    44.253825261852789 -18.769412812935450
    44.048416161483097 -18.381237814833021
    44.017390304153977 -18.081276075740377
    44.000984369468853 -17.861344565892828
    43.924889227540213 -17.552129101602546
    43.939133559415168 -17.501437246352943
    44.152405718021285 -17.124283107409664
    44.339700273404667 -16.854519712689463
    44.429478097211089 -16.659445671892001
    44.464576403540242 -16.193424052452059
    44.692595702572078 -16.172241632787330];

africa0_coord=[
    32.478361000000000 29.940417000000000
    32.625388999999998 28.976638999999999
    33.586221999999999 27.847500000000000
    35.138750000000002 24.513332999999999
    35.793749999999996 23.909222000000000
    35.476222000000000 23.934221999999998
    35.667888999999995 22.946610999999997
    36.895416999999995 22.063333000000000
    37.318750000000001 21.064999999999998
    37.111193999999998 21.214167000000000
    37.412056000000000 18.875000000000000
    38.557055999999996 18.060110999999999
    39.712027999999997 15.094999999999999
    39.875000000000000 15.507888999999999
    40.209972000000000 14.948722000000000
    41.189971999999997 14.623721999999999
    43.327027999999999 12.475778000000000
    43.376193999999998 11.988360999999999
    42.527471999999996 11.522888999999999
    43.156610999999998 11.612943999999999
    44.575806000000000 10.375389000000000
    45.805833000000000 10.867889000000000
    46.440805999999995 10.683750000000000
    47.412499999999994 11.179582999999999
    49.426666999999995 11.337888999999999
    50.797610999999996 11.986110999999999
    51.291249999999998 11.833306000000000
    51.022055999999999 10.420805999999999
    51.415416999999998 10.443332999999999
    50.902055999999995 10.314943999999999
    50.830416999999997 9.419167000000000
    47.965443999999998 4.477500000000000
    46.011666999999996 2.429528000000000
    43.559166999999995 0.714583000000000
    41.314971999999997 -1.961250000000000
    40.774555999999997 -1.934167000000000
    40.987916999999996 -2.264194000000000
    40.186278000000001 -2.735833000000000
    39.201943999999997 -4.680861000000000
    38.772888999999999 -6.052528000000000
    39.551249999999996 -7.017500000000000
    39.217110999999996 -7.850833000000000
    39.447055999999996 -7.813361000000000
    39.255389000000001 -8.292527999999999
    39.704443999999995 -10.031832999999999
    40.642055999999997 -10.684277999999999
    40.332082999999997 -11.313333000000000
    40.642083000000000 -12.746694000000000
    40.394556000000001 -12.928388999999999
    40.842110999999996 -14.832528000000000
    39.856639000000001 -16.438806000000000
    36.854971999999997 -17.878750000000000
    36.312500000000000 -18.877110999999999
    34.884999999999998 -19.852944000000001
    34.631222000000001 -19.650832999999999
    34.661221999999995 -20.555056000000000
    35.546222000000000 -22.171721999999999
    35.502055999999996 -24.107527999999999
    32.498666999999998 -25.960055999999998
    32.834055999999997 -26.292916999999999
    32.964388999999997 -26.084111000000000
    32.385388999999996 -28.549194000000000
    31.274555999999997 -29.454193999999998
    30.009166999999998 -31.297110999999997
    27.011638999999999 -33.565416999999997
    24.839971999999999 -34.206277999999998
    22.555833000000000 -33.978749999999998
    20.019971999999999 -34.830444000000000
    18.790806000000000 -34.080444000000000
    18.399943999999998 -34.302110999999996
    17.847888999999999 -32.820028000000001
    18.327082999999998 -32.530028000000001
    18.200443999999997 -31.682499999999997
    14.946249999999999 -26.317527999999999
    14.502027999999999 -22.542610999999997
    11.799500000000000 -17.979194000000000
    11.744555999999999 -15.835832999999999
    12.532888999999999 -13.418360999999999
    13.657055999999999 -12.235028000000000
    13.862833000000000 -11.008360999999999
    12.994582999999999 -9.087527999999999
    13.392111000000000 -8.388332999999999
    12.272917000000000 -6.140000000000000
    12.672943999999999 -6.010860999999999
    12.167916999999999 -5.727556000000000
    11.727860999999999 -4.450861000000000
    9.619583000000000 -2.408361000000000
    8.700388999999999 -0.668361000000000
    9.130417000000000 -0.714194000000000
    9.351638999999999 0.357917000000000
    9.501666999999999 0.050389000000000
    9.987083000000000 0.172500000000000
    9.309583000000000 0.512500000000000
    9.663750000000000 0.452500000000000
    9.695971999999999 1.078750000000000
    9.343722000000000 1.175833000000000
    9.975417000000000 3.079111000000000
    9.542916999999999 3.810833000000000
    9.833693999999999 3.809167000000000
    9.727943999999999 4.097500000000000
    9.214167000000000 3.941194000000000
    8.830389000000000 4.746667000000000
    8.554971999999999 4.482056000000000
    8.235778000000000 4.928750000000000
    8.278278000000000 4.536222000000000
    7.340000000000000 4.443722000000000
    7.075806000000000 4.737916999999999
    7.005777999999999 4.372917000000000
    6.779139000000000 4.667860999999999
    6.867471999999999 4.359583000000000
    6.723332999999999 4.612916999999999
    6.705028000000000 4.338722000000000
    5.983333000000000 4.310417000000000
    5.362083000000000 5.163306000000000
    5.680444000000000 5.540806000000000
    5.361667000000000 5.390416999999999
    5.176250000000000 5.568306000000000
    5.542860999999999 5.619193999999999
    5.080388999999999 5.718332999999999
    5.367083000000000 5.970000000000000
    5.064166999999999 5.764527999999999
    4.430833000000000 6.350389000000000
    3.401250000000000 6.395805999999999
    3.821639000000000 6.621194000000000
    1.343361000000000 6.161167000000000
    -1.981694000000005 4.754611000000000
    -4.016999999999996 5.258389000000000
    -6.124167000000000 4.997027999999999
    -7.480806000000030 4.360361000000000
    -8.717500000000030 4.807083000000000
    -11.230750000000000 6.789499999999999
    -12.521277999999995 7.384111000000000
    -12.456250000000011 7.766639000000000
    -13.293778000000032 8.421638999999999
    -12.837083000000007 8.560832999999999
    -13.174611000000027 8.531666999999999
    -13.303944000000001 9.334861000000000
    -13.723749999999995 9.506639000000000
    -13.537056000000007 9.796611000000000
    -14.455833000000041 10.208722000000000
    -14.697083000000021 11.076639000000000
    -15.078389000000016 10.918721999999999
    -15.004638999999997 11.183306000000000
    -15.235000000000014 10.997916999999999
    -15.504611000000011 11.328332999999999
    -15.033749999999998 11.643277999999999
    -15.509583000000021 11.789166999999999
    -14.939166999999998 11.988721999999999
    -15.955833000000041 11.734556000000000
    -15.853778000000034 11.992471999999999
    -16.352082999999993 12.096667000000000
    -16.014194000000032 12.345388999999999
    -16.776250000000005 12.406639000000000
    -16.220472000000029 12.609971999999999
    -16.763778000000002 12.571610999999999
    -16.717943999999989 13.459999999999999
    -16.355028000000004 13.249582999999999
    -15.574611000000004 13.515806000000000
    -16.506666999999993 13.361250000000000
    -16.513778000000002 13.879166999999999
    -17.533778000000041 14.741667000000000
    -16.577056000000027 15.697471999999999
    -16.043749999999989 17.747499999999999
    -16.547917000000041 19.379999999999999
    -16.197917000000018 20.224971999999998
    -16.912528000000009 21.161193999999998
    -17.107222000000036 20.837527999999999
    -16.983000000000004 21.798138999999999
    -15.704583000000014 23.972500000000000
    -15.915416999999991 23.803305999999999
    -14.903777999999988 24.685860999999999
    -14.465417000000002 26.196666999999998
    -13.651250000000005 26.659193999999999
    -12.915861000000007 27.956249999999997
    -11.525000000000034 28.298749999999998
    -10.217943999999989 29.309138999999998
    -9.662944000000039 30.071638999999998
    -9.846222000000012 31.397499999999997
    -9.287111000000039 32.545833000000002
    -6.812500000000000 34.025388999999997
    -5.934639000000004 35.790777999999996
    -5.470861000000014 35.917110999999998
    -4.688389000000029 35.206277999999998
    -2.979139000000032 35.441972000000000
    -2.830083000000002 35.108666999999997
    -1.957500000000039 35.075389000000001
    -0.937528000000043 35.727916999999998
    0.000000000000000 35.842889000000000
    1.084139000000000 36.493749999999999
    3.895861000000000 36.925388999999996
    5.344139000000000 36.641193999999999
    6.415806000000000 37.086250000000000
    7.899972000000000 36.842917000000000
    9.744999999999999 37.349556000000000
    10.431611000000000 36.712888999999997
    11.039138999999999 37.089582999999998
    10.467832999999999 36.135777999999995
    11.162917000000000 35.215832999999996
    10.017861000000000 34.187528000000000
    10.214138999999999 33.797083000000001
    12.299139000000000 32.837888999999997
    13.388306000000000 32.897860999999999
    15.205027999999999 32.382055999999999
    15.741667000000000 31.397888999999999
    17.853306000000000 30.925388999999999
    18.989971999999998 30.274583000000000
    20.057082999999999 30.850860999999998
    20.099999999999998 32.182055999999996
    21.708306000000000 32.941221999999996
    23.107471999999998 32.635388999999996
    23.319944000000000 32.152889000000002
    24.979832999999999 31.971055999999997
    25.200832999999999 31.523638999999999
    26.982471999999998 31.442055999999997
    29.064943999999997 30.822056000000000
    31.028305999999997 31.600389000000000
    31.919166999999998 31.531222000000000
    32.478361000000000 29.940417000000000];

if flag
    madagascar_pshape=polyshape(madagascar_coord);
    africa0_pshape=polyshape(africa0_coord);
    africa_pshape=union(madagascar_pshape,africa0_pshape);
else
    africa_pshape=polyshape(africa0_coord);
end