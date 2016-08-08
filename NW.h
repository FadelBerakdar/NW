/*
 * =====================================================================================
 *
 *       Filename:  alignment.h
 *    Description:  Alignment Module 
 *	      URL:  http://www3.cs.stonybrook.edu/~rp/class/549f14/lectures/CSE549-Lec04.pdf
 *        Version:  1.0
 *        Created:  08/02/2016 14:37:30
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Fadel Berakdar
 *   Organization:  Novocraft
 *
 * =====================================================================================
 */

#ifndef NW_H
#define NW_H


class Alignment {
    
    private: 
	int  score, match_mismatch, i, j, k, m, n, indexI, indexJ, indexK, sigma, epsilon, negative;
	const char *read, *target;
    	int  **V, **G, **E, **F;			
	char **tracingBack;
    	char figure;
	
	int  max(int diagonal,  int vertical, int horizontal, char *figure) ;
    	int  substitutionValue(char a, char b); 	

    public:
	Alignment(const char *_target, const char *_read, int _sigma, int _epsilon);
	~Alignment();
	void global() ;

	static int  substitutionMatrix[4][4]; 

};

#endif



