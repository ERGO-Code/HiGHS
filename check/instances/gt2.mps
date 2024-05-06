*NAME:          gt2
*ROWS:          29
*COLUMNS:       188
*INTEGER:       188
*NONZERO:       376
*BEST SOLN:     21166.000
*LP SOLN:       13460.233074
*SOURCE:
*		Laurence A. Wolsey (University of Louvain)
*               Sebastian Ceria (Columbia University)
*APPLICATION:   Truck routing
*COMMENTS:      24 of the integer variables are binary
*
*
NAME          GT2 
ROWS
 N  COST....
 G  dem...01
 G  dem...02
 G  dem...03
 G  dem...04
 G  dem...05
 G  dem...06
 G  dem...07
 G  dem...08
 G  dem...09
 G  dem...10
 G  dem...11
 G  dem...12
 L  avail.01
 L  avail.02
 L  avail.03
 L  avail.04
 L  avail.05
 L  avail.06
 L  avail.07
 L  avail.08
 L  avail.09
 L  avail.10
 L  avail.11
 L  avail.12
 L  avail.13
 L  avail.14
 L  avail.15
 L  avail.16
 L  avail.17
COLUMNS
    MARK0000  'MARKER'                 'INTORG'
    x...0101  COST....          3595   dem...01            50
    x...0101  avail.01             1
    x...0201  COST....          3595   dem...02            50
    x...0201  avail.01             1
    x...0301  COST....          1833   dem...03            50
    x...0301  avail.01             1
    x...0401  COST....          3595   dem...04            50
    x...0401  avail.01             1
    x...0701  COST....          3595   dem...07            50
    x...0701  avail.01             1
    x...0801  COST....          3595   dem...08            50
    x...0801  avail.01             1
    x...0901  COST....          1833   dem...09            50
    x...0901  avail.01             1
    x...1001  COST....          3595   dem...10            50
    x...1001  avail.01             1
    x...0102  COST....          3100   dem...01            45
    x...0102  avail.02             1
    x...0202  COST....          3100   dem...02            45
    x...0202  avail.02             1
    x...0302  COST....          1257   dem...03            45
    x...0302  avail.02             1
    x...0402  COST....          3100   dem...04            45
    x...0402  avail.02             1
    x...0502  COST....          1355   dem...05            28
    x...0502  avail.02             1
    x...0602  COST....          7797   dem...06            28
    x...0602  avail.02             1
    x...0702  COST....          3100   dem...07            45
    x...0702  avail.02             1
    x...0802  COST....          3100   dem...08            45
    x...0802  avail.02             1
    x...0902  COST....          1257   dem...09            45
    x...0902  avail.02             1
    x...1002  COST....          3100   dem...10            45
    x...1002  avail.02             1
    x...1102  COST....          1355   dem...11            28
    x...1102  avail.02             1
    x...1202  COST....          7797   dem...12            28
    x...1202  avail.02             1
    x...0103  COST....          2685   dem...01            40
    x...0103  avail.03             1
    x...0203  COST....          2685   dem...02            40
    x...0203  avail.03             1
    x...0303  COST....          1257   dem...03            40
    x...0303  avail.03             1
    x...0403  COST....          2685   dem...04            40
    x...0403  avail.03             1
    x...0503  COST....          1176   dem...05            45
    x...0503  avail.03             1
    x...0603  COST....          6214   dem...06            45
    x...0603  avail.03             1
    x...0703  COST....          2685   dem...07            40
    x...0703  avail.03             1
    x...0803  COST....          2685   dem...08            40
    x...0803  avail.03             1
    x...0903  COST....          1257   dem...09            40
    x...0903  avail.03             1
    x...1003  COST....          2685   dem...10            40
    x...1003  avail.03             1
    x...1103  COST....          1176   dem...11            45
    x...1103  avail.03             1
    x...1203  COST....          6214   dem...12            45
    x...1203  avail.03             1
    x...0104  COST....          2940   dem...01            30
    x...0104  avail.04             1
    x...0204  COST....          2940   dem...02            30
    x...0204  avail.04             1
    x...0304  COST....          1980   dem...03            30
    x...0304  avail.04             1
    x...0404  COST....          2940   dem...04            30
    x...0404  avail.04             1
    x...0504  COST....          1980   dem...05          16.5
    x...0504  avail.04             1
    x...0604  COST....          7145   dem...06          16.5
    x...0604  avail.04             1
    x...0704  COST....          2940   dem...07            30
    x...0704  avail.04             1
    x...0804  COST....          2940   dem...08            30
    x...0804  avail.04             1
    x...0904  COST....          1980   dem...09            30
    x...0904  avail.04             1
    x...1004  COST....          2940   dem...10            30
    x...1004  avail.04             1
    x...1104  COST....          1980   dem...11          16.5
    x...1104  avail.04             1
    x...1204  COST....          7145   dem...12          16.5
    x...1204  avail.04             1
    x...0105  COST....          3100   dem...01            32
    x...0105  avail.05             1
    x...0205  COST....          3100   dem...02            32
    x...0205  avail.05             1
    x...0305  COST....          1816   dem...03            32
    x...0305  avail.05             1
    x...0405  COST....          3100   dem...04            32
    x...0405  avail.05             1
    x...0505  COST....          1926   dem...05            45
    x...0505  avail.05             1
    x...0605  COST....          7011   dem...06            45
    x...0605  avail.05             1
    x...0705  COST....          3100   dem...07            32
    x...0705  avail.05             1
    x...0805  COST....          3100   dem...08            32
    x...0805  avail.05             1
    x...0905  COST....          1816   dem...09            32
    x...0905  avail.05             1
    x...1005  COST....          3100   dem...10            32
    x...1005  avail.05             1
    x...1105  COST....          1926   dem...11            45
    x...1105  avail.05             1
    x...1205  COST....          7011   dem...12            45
    x...1205  avail.05             1
    x...0106  COST....          3100   dem...01            50
    x...0106  avail.06             1
    x...0206  COST....          3100   dem...02            50
    x...0206  avail.06             1
    x...0306  COST....          1816   dem...03            50
    x...0306  avail.06             1
    x...0406  COST....          3100   dem...04            52
    x...0406  avail.06             1
    x...0506  COST....          1983   dem...05            52
    x...0506  avail.06             1
    x...0606  COST....          7797   dem...06            52
    x...0606  avail.06             1
    x...0706  COST....          3100   dem...07            50
    x...0706  avail.06             1
    x...0806  COST....          3100   dem...08            50
    x...0806  avail.06             1
    x...0906  COST....          1816   dem...09            50
    x...0906  avail.06             1
    x...1006  COST....          3100   dem...10            52
    x...1006  avail.06             1
    x...1106  COST....          1983   dem...11            52
    x...1106  avail.06             1
    x...1206  COST....          7797   dem...12            52
    x...1206  avail.06             1
    x...0507  COST....          1926   dem...05            45
    x...0507  avail.07             1
    x...0607  COST....          7011   dem...06            45
    x...0607  avail.07             1
    x...1107  COST....          1926   dem...11            45
    x...1107  avail.07             1
    x...1207  COST....          7011   dem...12            45
    x...1207  avail.07             1
    x...0108  COST....          3100   dem...01            45
    x...0108  avail.08             1
    x...0208  COST....          3100   dem...02            45
    x...0208  avail.08             1
    x...0308  COST....          1257   dem...03            45
    x...0308  avail.08             1
    x...0408  COST....          3100   dem...04            44
    x...0408  avail.08             1
    x...0708  COST....          3100   dem...07            45
    x...0708  avail.08             1
    x...0808  COST....          3100   dem...08            45
    x...0808  avail.08             1
    x...0908  COST....          1257   dem...09            45
    x...0908  avail.08             1
    x...1008  COST....          3100   dem...10            44
    x...1008  avail.08             1
    x...0109  COST....          2448   dem...01          2534
    x...0109  avail.09             1
    x...0209  COST....          2448   dem...02          2534
    x...0209  avail.09             1
    x...0309  COST....          1652   dem...03          2534
    x...0309  avail.09             1
    x...0409  COST....          2448   dem...04          2534
    x...0409  avail.09             1
    x...0509  COST....          1652   dem...05          2519
    x...0509  avail.09             1
    x...0609  COST....          5954   dem...06          2519
    x...0609  avail.09             1
    x...0709  COST....          2448   dem...07          2534
    x...0709  avail.09             1
    x...0809  COST....          2448   dem...08          2534
    x...0809  avail.09             1
    x...0909  COST....          1652   dem...09          2534
    x...0909  avail.09             1
    x...1009  COST....          2448   dem...10          2534
    x...1009  avail.09             1
    x...1109  COST....          1652   dem...11          2519
    x...1109  avail.09             1
    x...1209  COST....          5954   dem...12          2519
    x...1209  avail.09             1
    x...0110  dem...01            22   avail.10             1
    x...0111  dem...01            56   avail.11             1
    x...0112  dem...01            30   avail.12             1
    x...0113  dem...01            30   avail.13             1
    x...0114  dem...01            23   avail.14             1
    x...0115  dem...01            42   avail.15             1
    x...0116  dem...01            56   avail.16             1
    x...0117  dem...01            28   avail.17             1
    x...0210  dem...02            22   avail.10             1
    x...0211  dem...02            56   avail.11             1
    x...0212  dem...02            30   avail.12             1
    x...0213  dem...02            30   avail.13             1
    x...0214  dem...02            23   avail.14             1
    x...0215  dem...02            42   avail.15             1
    x...0216  dem...02            56   avail.16             1
    x...0217  dem...02            28   avail.17             1
    x...0310  dem...03            22   avail.10             1
    x...0311  dem...03            56   avail.11             1
    x...0312  dem...03            30   avail.12             1
    x...0313  dem...03            30   avail.13             1
    x...0314  dem...03            23   avail.14             1
    x...0315  dem...03            42   avail.15             1
    x...0316  dem...03            56   avail.16             1
    x...0317  dem...03            28   avail.17             1
    x...0410  dem...04            22   avail.10             1
    x...0411  dem...04            56   avail.11             1
    x...0412  dem...04            30   avail.12             1
    x...0413  dem...04            30   avail.13             1
    x...0414  dem...04            23   avail.14             1
    x...0415  dem...04            42   avail.15             1
    x...0416  dem...04            56   avail.16             1
    x...0417  dem...04            28   avail.17             1
    x...0510  dem...05            25   avail.10             1
    x...0511  dem...05            34   avail.11             1
    x...0512  dem...05            25   avail.12             1
    x...0513  dem...05            19   avail.13             1
    x...0514  dem...05            56   avail.14             1
    x...0515  dem...05            34   avail.15             1
    x...0516  dem...05            34   avail.16             1
    x...0517  dem...05            25   avail.17             1
    x...0610  dem...06            25   avail.10             1
    x...0611  dem...06            34   avail.11             1
    x...0612  dem...06            25   avail.12             1
    x...0613  dem...06            19   avail.13             1
    x...0614  dem...06            56   avail.14             1
    x...0615  dem...06            34   avail.15             1
    x...0616  dem...06            34   avail.16             1
    x...0617  dem...06            25   avail.17             1
    x...0710  dem...07            22   avail.10             1
    x...0711  dem...07            56   avail.11             1
    x...0712  dem...07            30   avail.12             1
    x...0713  dem...07            30   avail.13             1
    x...0714  dem...07            23   avail.14             1
    x...0715  dem...07            42   avail.15             1
    x...0716  dem...07            56   avail.16             1
    x...0717  dem...07            28   avail.17             1
    x...0810  dem...08            22   avail.10             1
    x...0811  dem...08            56   avail.11             1
    x...0812  dem...08            30   avail.12             1
    x...0813  dem...08            30   avail.13             1
    x...0814  dem...08            23   avail.14             1
    x...0815  dem...08            42   avail.15             1
    x...0816  dem...08            56   avail.16             1
    x...0817  dem...08            28   avail.17             1
    x...0910  dem...09            22   avail.10             1
    x...0911  dem...09            56   avail.11             1
    x...0912  dem...09            30   avail.12             1
    x...0913  dem...09            30   avail.13             1
    x...0914  dem...09            23   avail.14             1
    x...0915  dem...09            42   avail.15             1
    x...0916  dem...09            56   avail.16             1
    x...0917  dem...09            28   avail.17             1
    x...1010  dem...10            22   avail.10             1
    x...1011  dem...10            56   avail.11             1
    x...1012  dem...10            30   avail.12             1
    x...1013  dem...10            30   avail.13             1
    x...1014  dem...10            23   avail.14             1
    x...1015  dem...10            42   avail.15             1
    x...1016  dem...10            56   avail.16             1
    x...1017  dem...10            28   avail.17             1
    x...1110  dem...11            25   avail.10             1
    x...1111  dem...11            34   avail.11             1
    x...1112  dem...11            25   avail.12             1
    x...1113  dem...11            19   avail.13             1
    x...1114  dem...11            56   avail.14             1
    x...1115  dem...11            34   avail.15             1
    x...1116  dem...11            34   avail.16             1
    x...1117  dem...11            25   avail.17             1
    x...1210  dem...12            25   avail.10             1
    x...1211  dem...12            34   avail.11             1
    x...1212  dem...12            25   avail.12             1
    x...1213  dem...12            19   avail.13             1
    x...1214  dem...12            56   avail.14             1
    x...1215  dem...12            34   avail.15             1
    x...1216  dem...12            34   avail.16             1
    x...1217  dem...12            25   avail.17             1
    MARK0001  'MARKER'                 'INTEND'
RHS
    rhs       dem...01           200   dem...02            80
    rhs       dem...03           300   dem...04           250
    rhs       dem...05            61   dem...06          6064
    rhs       dem...07            94   dem...08            67
    rhs       dem...09           450   dem...10           231
    rhs       dem...11            76   avail.01             9
    rhs       avail.02            15   avail.03            14
    rhs       avail.04             1   avail.05             6
    rhs       avail.06             8   avail.07             5
    rhs       avail.08            12   avail.09             5
    rhs       avail.10             6   avail.11             4
    rhs       avail.12             9   avail.13             5
    rhs       avail.14             4   avail.15             2
    rhs       avail.16             2   avail.17             1
BOUNDS
 UP bnd       x...0101             9
 UP bnd       x...0201             9
 UP bnd       x...0301             9
 UP bnd       x...0401             9
 UP bnd       x...0701             9
 UP bnd       x...0801             9
 UP bnd       x...0901             9
 UP bnd       x...1001             9
 UP bnd       x...0102            15
 UP bnd       x...0202            15
 UP bnd       x...0302            15
 UP bnd       x...0402            15
 UP bnd       x...0502            15
 UP bnd       x...0602            15
 UP bnd       x...0702            15
 UP bnd       x...0802            15
 UP bnd       x...0902            15
 UP bnd       x...1002            15
 UP bnd       x...1102            15
 UP bnd       x...1202            15
 UP bnd       x...0103            14
 UP bnd       x...0203            14
 UP bnd       x...0303            14
 UP bnd       x...0403            14
 UP bnd       x...0503            14
 UP bnd       x...0603            14
 UP bnd       x...0703            14
 UP bnd       x...0803            14
 UP bnd       x...0903            14
 UP bnd       x...1003            14
 UP bnd       x...1103            14
 UP bnd       x...1203            14
 UP bnd       x...0104             1
 UP bnd       x...0204             1
 UP bnd       x...0304             1
 UP bnd       x...0404             1
 UP bnd       x...0504             1
 UP bnd       x...0604             1
 UP bnd       x...0704             1
 UP bnd       x...0804             1
 UP bnd       x...0904             1
 UP bnd       x...1004             1
 UP bnd       x...1104             1
 UP bnd       x...1204             1
 UP bnd       x...0105             6
 UP bnd       x...0205             6
 UP bnd       x...0305             6
 UP bnd       x...0405             6
 UP bnd       x...0505             6
 UP bnd       x...0605             6
 UP bnd       x...0705             6
 UP bnd       x...0805             6
 UP bnd       x...0905             6
 UP bnd       x...1005             6
 UP bnd       x...1105             6
 UP bnd       x...1205             6
 UP bnd       x...0106             8
 UP bnd       x...0206             8
 UP bnd       x...0306             8
 UP bnd       x...0406             8
 UP bnd       x...0506             8
 UP bnd       x...0606             8
 UP bnd       x...0706             8
 UP bnd       x...0806             8
 UP bnd       x...0906             8
 UP bnd       x...1006             8
 UP bnd       x...1106             8
 UP bnd       x...1206             8
 UP bnd       x...0507             5
 UP bnd       x...0607             5
 UP bnd       x...1107             5
 UP bnd       x...1207             5
 UP bnd       x...0108            12
 UP bnd       x...0208            12
 UP bnd       x...0308            12
 UP bnd       x...0408            12
 UP bnd       x...0708            12
 UP bnd       x...0808            12
 UP bnd       x...0908            12
 UP bnd       x...1008            12
 UP bnd       x...0109             5
 UP bnd       x...0209             5
 UP bnd       x...0309             5
 UP bnd       x...0409             5
 UP bnd       x...0509             5
 UP bnd       x...0609             5
 UP bnd       x...0709             5
 UP bnd       x...0809             5
 UP bnd       x...0909             5
 UP bnd       x...1009             5
 UP bnd       x...1109             5
 UP bnd       x...1209             5
 UP bnd       x...0110             6
 UP bnd       x...0111             4
 UP bnd       x...0112             9
 UP bnd       x...0113             5
 UP bnd       x...0114             4
 UP bnd       x...0115             2
 UP bnd       x...0116             2
 UP bnd       x...0117             1
 UP bnd       x...0210             6
 UP bnd       x...0211             4
 UP bnd       x...0212             9
 UP bnd       x...0213             5
 UP bnd       x...0214             4
 UP bnd       x...0215             2
 UP bnd       x...0216             2
 UP bnd       x...0217             1
 UP bnd       x...0310             6
 UP bnd       x...0311             4
 UP bnd       x...0312             9
 UP bnd       x...0313             5
 UP bnd       x...0314             4
 UP bnd       x...0315             2
 UP bnd       x...0316             2
 UP bnd       x...0317             1
 UP bnd       x...0410             6
 UP bnd       x...0411             4
 UP bnd       x...0412             9
 UP bnd       x...0413             5
 UP bnd       x...0414             4
 UP bnd       x...0415             2
 UP bnd       x...0416             2
 UP bnd       x...0417             1
 UP bnd       x...0510             6
 UP bnd       x...0511             4
 UP bnd       x...0512             9
 UP bnd       x...0513             5
 UP bnd       x...0514             4
 UP bnd       x...0515             2
 UP bnd       x...0516             2
 UP bnd       x...0517             1
 UP bnd       x...0610             6
 UP bnd       x...0611             4
 UP bnd       x...0612             9
 UP bnd       x...0613             5
 UP bnd       x...0614             4
 UP bnd       x...0615             2
 UP bnd       x...0616             2
 UP bnd       x...0617             1
 UP bnd       x...0710             6
 UP bnd       x...0711             4
 UP bnd       x...0712             9
 UP bnd       x...0713             5
 UP bnd       x...0714             4
 UP bnd       x...0715             2
 UP bnd       x...0716             2
 UP bnd       x...0717             1
 UP bnd       x...0810             6
 UP bnd       x...0811             4
 UP bnd       x...0812             9
 UP bnd       x...0813             5
 UP bnd       x...0814             4
 UP bnd       x...0815             2
 UP bnd       x...0816             2
 UP bnd       x...0817             1
 UP bnd       x...0910             6
 UP bnd       x...0911             4
 UP bnd       x...0912             9
 UP bnd       x...0913             5
 UP bnd       x...0914             4
 UP bnd       x...0915             2
 UP bnd       x...0916             2
 UP bnd       x...0917             1
 UP bnd       x...1010             6
 UP bnd       x...1011             4
 UP bnd       x...1012             9
 UP bnd       x...1013             5
 UP bnd       x...1014             4
 UP bnd       x...1015             2
 UP bnd       x...1016             2
 UP bnd       x...1017             1
 UP bnd       x...1110             6
 UP bnd       x...1111             4
 UP bnd       x...1112             9
 UP bnd       x...1113             5
 UP bnd       x...1114             4
 UP bnd       x...1115             2
 UP bnd       x...1116             2
 UP bnd       x...1117             1
 UP bnd       x...1210             6
 UP bnd       x...1211             4
 UP bnd       x...1212             9
 UP bnd       x...1213             5
 UP bnd       x...1214             4
 UP bnd       x...1215             2
 UP bnd       x...1216             2
 UP bnd       x...1217             1
ENDATA
