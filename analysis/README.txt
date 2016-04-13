Progam overview as implemented in runFullAnalysis.sh

():progamms in ( ) are optional drawing functions
{}:special requirments, IR=inclusive region
*:comments

*************************************************************************  ************************  *******************************
*                 BABY TREES SKIMMED                                    *  * BABY TREES QCD SKIM  *  * BABY TREES QCD MONOJET SKIM *
*************************************************************************  ************************  *******************************
   |              |                  |                  |           |			|			|
   |              |                  |                  |           V			|			|
   |              |                  |	                |     zllControlRegion		|			|
   |              |                  |                  |    (drawZllControlRegion)	|			|
   |              |                  |                  |	        |{IR}		|			|
   |              V                  |                  V	        |		|			|
   |         regionEventYields	     |          makeZinvGammaTemplates  |   		|			|
   |             |    |      |	     |	       	(drawGammaTemplates)    |		|			|
   |             |    |      |	     |	       	        |{IR}	        |		|			|
   |             |    |      |	     |                  |	        |		|			|
   |             |    |      |	     |                  |	        |		|			|
   |             |    |      |	     V                  |	        |		|			V
   |             |    |      |  gammaControlRegion------|-------\ 	|		|	     qcdControlRegion (monojet)	                 
   |             |    |      |  (drawGammaControlRegion)|       | 	|		|		 	|{qcd IR}
   |             |    |      |       |       |          |       | 	|		|		 	|
   |             |    |      |       |       |          |       | 	|		|		 	|
   |             |    |      |       |       |          |       | 	|		|		 	|
   |             |    |      |       |       |          |       | 	|		V		 	|
   |             |    |      |       |       |          |       |	|	qcdControlRegion 		|
   | 	         |    |      |       |       V          |       |	|	        |{qcd IR}		|	     
   | 	     	 |    |      |       | binIsoAlongAxes  |       |	|	        |			|
   | 		 |    |      |       |       |{IR}      |       |	|		|			|
   | 		 |    |      |       |       |          |       |	|		|			|
   | 		 |    |      |       |       |          |       |	|		|			|
   | 		 |    |      |       |       |          |       |	|		|			|
   | 		 |    |      |       |       |          |       |	|		|			V
   V 		 |    |      |       |       V	        V       |	|		|	        drawQCDForMonoJet
llControlRegion  |    |      |       |	   fitPurityGamma       |	|		|			|
    	|	 |    |      |       |    (drawPurityFits)      |	|		V			|
    	|	 |    |      |       |       |{IR}      |       |	|     computeQCDFromDeltaPhi		|
    	|	 |    |      |       |       |          |	|	|     (drawQCDFromDeltaPhi)		|
    	|	 |    |      |       |       |          |       |	|		|			|
    	|	 |    |      |       |       |          |       |	|		|			|
    	|	 |    |      |       |       |          |       |	|		|			|
    	|	 |    |      |       |       |          |       |	|		|			|
        V     	 V    |      V       V       V	        V	V	V		|			|
    computeLostLepton |     computeZinvFromGamma      computeZllGammaRatio		|			|
        |             |       	|		         	|{IR}			|			|
        |             |       	|		             	|			|			|
        |             |       	|		        	|			|			|
        |             |       	|		        	|			|			|
	|	      |	 	\_______________________________/			|			|
	|	      \_________________________|_______________________________________/			|
	\_______________________________________|_______________________________________________________________/				
					        |
					  createDatacards
						|								     
						|		
        					|	
        					V
	       	  ****************************************************************     
		  *                       DATACARDS                              *     
		  ****************************************************************    

