//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

const char * welcomeScreen = R"(
                                       MMMMM
                                      M8NO  +M
                                    NO       M         ,MM
                                    MM     ?MN,   D8MMMO  :NMM
                               :MZ            ZM  MZ            NM,  M
                              M,          :M   MMN            DMMN M  ~
                             M        8MMM M  M,    7O: :M, M+     N:
                            M  =NMZ M     ,M8 M    ,O     MM   :    M
                            D  DM   M    $8        I= MM  $M   7    M
                           Z~ M M  O 8   8          M     N MM    DM
                           8, O    MM    M           :MM8,     ,, MMD
                           O:       M    87MMO ,  MM~                ,M,
                           I=O  M   M   M        M                     M
                           :M   M   N   M N,     M                    M~
                           :M  ?~  N    8I+  M    MMMMMMMMMMMMMMMOM? 7,
                            M  ++  M   8MZ7  M     M                MM
                            M   D  M  N   N   M,    ND          DM8,
                            M   M  M M~M               DMNM88= 8,
                            M    8  ?M?       M                 M~
                            M      N~8       ?$    :             M
                            NN      M     M  M     M             8:
                            D M    M   DDM   M     M               M,
                            D     ,D  D?M    N     M                 ND
                            M     M  ,O8$   ,7     M    M    M         M~
                            NN  8IO  M M    I=     =    M    M          ,M
                                 M  M                   M    M I          M
                                ~M ?8                 M M    M +M   MM     N
        N+I8MMM~                MM M ,      M         M D    M   M   MON   ~$
        M       NM              M8   M      O     D   M Z+  M:    M,  M M   M
         ,MMMD    88            8M   N     M      M         M      ~M  M M+  M
               7M, MI          +ZM  M      M      M        I$ $NI    M ,O M8 ?O
                 ?N M          I:   M      M       8       M     ,    M M     M
                  M +?         8,          N       M       N    :MMMMM7D M    M
                   D M         8,                  N,     ,$   M       M 8    8:
                   N,D,        ?:    =MMMMMZ               D  M         M  M  8ZMM
                    M ND       :Z=MO,                      N    MM8ZN8  M  M  M+ ,:~MMM 8~
                     M  DMMMMMNZ             D    ~MM        :MM        M   M MM    MNMM M
                     $7           DMMZZNOMMNMZ       M     =M       MZ  :7 MM   MM$ M8   N$ M
                       MO      MMZ        Z                +   ,MMD   +  D M~8  N  M       M
                          ?I=    D         8                        ,NMMMMO MM M 8  M     ID
                                  8        D                     MN,     M  ~,  MM       8N
     ,MNZ,  ,8MM         :~ZMMM,  :8       D                   8         ,M      N
   88           :M   M+         M  ,D       D                $            M      M
  M:       M      M ,7           D  8N$     8~             :D             M      M ?M,
  M    7   M      M M   +MM?     M ,8  D8    D=          ,D              Z?     ,O
  $I~ ~?  Z8     M  M      M     M 7+  M  :,   D$      8Z  M? M          M$D+M  7M?
   Z:  MMN       N, M     M?     M D   M          =$~   M :+  D,      NMMMMM OMM  MM          MMM
   O:           :M $Z    D8      M M   M          M     M O:  N      M     ,O     M M        M   N
   $~          M:  M  D  DM      O N   M   MM     M    :N M   M     +7      M     M  M       D   D
   ~=  MNMMMM8    D: $     M    ?Z M   D   M M    M    M  M   M     M M M   7,    M   M     D~   7=
   ~I     M       N  O     :$   N  M   D   8 Z=   M   M~  M   M     Z   N    M    M   D,    M    Z
    M     M      M   D     $~   M  M   N   M  M   M   M   Z   M    M         D    M    N    M    M
    M     M      O   M     M    M  M   M  M     ,M   Z7  N~   N    M  NM     M    N    M    M    M
 M  M     M     M    M     D    M  M   MMM   M M     M   M    M    N M  N    M   7$     D   M    M
  MMM7  MN     O$   M     M    N   M          +      M   M    $   N +I  ~+   M   M  ,   M   M    M
                MMM+      M    M   M       M,        N   M   O=   M M   ~?   M   M M~NN  M  8    M
                           N  $$   7M      M+       D    8   M    M M   D,   M   M M  M  ,N 8    D
                            D8:     M      OM      MN    D   M  , M M   M    M   O $  +~  $MZ    D
                                     :M88DM  N +MM+      MMM$M  N~+~N   D,   M  :8  D  M         N
                                                                 MMM     M   M  Z=  M  M         D

//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

   )";

