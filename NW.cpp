
#include "NW.h"
#include <string>
#include <iostream>
#include <cstring>
#include <iomanip>



Alignment::Alignment(const char *_target,const char *_read, int _sigma, int _epsilon)
:sigma(_sigma), epsilon(_epsilon)
{
    read = _read;
    target = _target;
    m = strlen(read) + 1 ;
    n = strlen(target) + 1 ;
    negative = -10000;
    
    // Initilaizing of 3D dynamic Scoring matrix
    V = new int *[m];
    for (i = 0; i < m; i++){
        V[i] = new int[n] ;
    }
    
    G = new int *[m];
    for (i = 0; i < m; i++){
        G[i] = new int[n] ;
    }
    
    
    E = new int *[m];
    for (i = 0; i < m; i++){
        E[i] = new int[n] ;
    }
    
    
    F = new int *[m];
    for (i = 0; i < m; i++){
        F[i] = new int[n] ;
    }
    
    // std::cout << "Finished Initialization of scoring :... " << std::endl;
    
    // Initilaizing of 3D dynamic Back Tracking matrix
    tracingBack = new char *[m];
    for (i = 0; i < m; i++){
        tracingBack[i] = new char[n] ;
    }
}



Alignment::~Alignment() {
    // *** ********
    // deleting the scoring matrices:
    for (i =0; i < m ; i++) {
        delete [] V[i];
    }
    delete [] V;
    
    for (i =0; i < m ; i++) {
        delete [] G[i];
    }
    delete [] G;
    
    for (i =0; i < m ; i++) {
        delete [] E[i];
    }
    delete [] E;
    
    for (i =0; i < m ; i++) {
        delete [] F[i];
    }
    delete [] F;
    
    
    for (i =0; i < m ; i++) {
        delete [] tracingBack[i];
    }
    delete [] tracingBack;
    //std::cout << "Finished Deleting all matrics :... " << std::endl;
}


void Alignment::global() {
    std::string alignRef, alignTarget, dots;
    // filling the first row, first column in the scoring matrix and backtracking matrix
    for (i = 0; i < m ; i++ ) {
        for (j = 0 ; j < n ; j++ ) {
            
            if (i == 0 && j == 0) {
                V[0][0]  = 0;
                tracingBack[0][0] = 'n';
            }
            
            else if (i == 0) {
                V[0][j] = E[0][j] =  0;
                G[0][j] = F[0][j] = 0;
                tracingBack[0][j] = '-';
            }
            
            else if (j == 0) {
                V[i][0] = F[i][0] =  0;
                G[i][0] = E[i][0] = 0;
                tracingBack[i][0] = '|';
                
            }
        }
    }
    // std::cout << "Finished Filling up the first row & colum of both matrices:... " << std::endl;
    
    for (i = 1; i < m; i++)
    {
        for (j = 1; j< n; j++)
        {
            //match_mismatch = substitutionValue(read[i-1], target[j-1]);
            match_mismatch = read[i-1] == target[j-1] ?  5 : -4; // part of DNAFULL
            G[i][j] = V[i-1][j-1] + match_mismatch;				  //
            E[i][j] = std::max(F[i-1][j] - epsilon, V[i-1][j] - sigma - epsilon); // -
            F[i][j] = std::max(E[i][j-1] - epsilon, V[i][j-1] - sigma - epsilon); // |
            V[i][j] = max(G[i][j], E[i][j], F[i][j], &figure);
            tracingBack[i][j] = figure;
        }
        // std::cout << "Finished Filling up the rest rows & colums of both matrices:... " << std::endl;
    }
    

    // Global Aligment strings

    
    i = m - 1 , j = n - 1;
    score = negative;
    
    for ( j = 1; j < n; j++) {
       if  ( V[i][j] >= score ) {
                score = V[i][j];
                indexJ = j;
       }
    }
    
    i = m - 1 , j = indexJ;
    while ( tracingBack[i][j] != 'n' ){
        
        figure = tracingBack[i][j];
        
        switch(figure) {
                
            case '\\' :
                alignTarget = read[i-1] + alignTarget ;
                alignRef = target[j-1] + alignRef ;
                read[i-1] == target[j-1] ? dots = '|' + dots : dots = '.' + dots ;
                i--;
                j--;
                break;
                
                
            case '|' :
                alignTarget = read[i-1] + alignTarget ;
                alignRef = '-' + alignRef;
                dots = ' ' + dots;
                i--;
                break;
                
                
            case '-' :	    alignTarget = '-' + alignTarget ;
                alignRef = target[j-1] + alignRef;
                dots = ' ' + dots;
                j--;
        }
    }
    
    // adding up the rest of the target residues and - in read
    for ( j = indexJ  ; j < n -1 ; j++ ){
        alignRef += target[j];
        alignTarget += '-';
    }

    
    std::cout << "\n***  ";
    std::cout << "Global Alignment";
    std::cout << "  ***";
    std::cout << "\n" << alignRef << std::endl  << dots << std::endl << alignTarget << std::endl;
    std::cout << "Score is: " << score << std::endl;
}





int Alignment::substitutionValue(char a, char b) {
    
    // takes tow bases and returns integer equals to their value in substitution matrix
    // it can be extended to other characters like N,
    
    int x, y;
    
    switch(a) {
        case 'A': x = 0; break;
        case 'T': x = 1; break;
        case 'C': x = 2; break;
        case 'G': x = 3; break;
    }
    
    switch(b) {
        case 'A': y = 0; break;
        case 'T': y = 1; break;
        case 'C': y = 2; break;
        case 'G': y = 3; break;
    }
    
    
    return substitutionMatrix[x][y];
}


// part of DNAfull http://rosalind.info/glossary/dnafull/
//       A    T     C     G
int Alignment::substitutionMatrix[4][4] =  {  {  5,  -4,   -4,   -4},  //  A
    { -4,   5,   -4,   -4},  //  T
    { -4,  -4,    5,   -4},  //  C
    { -4,  -4,   -4,    5},  //  G
};



int Alignment::max(int diagonal,  int vertical, int horizontal, char *figure) {
    
    int  max = 0 ;
    
    if( diagonal > vertical && diagonal > horizontal )
    {
        max = diagonal ;
        *figure = '\\' ;
    }
    
    else if (vertical > horizontal)
    {
        max = vertical ;
        *figure = '|' ;
    }
    
    else     
    {
        max = horizontal ;
        *figure = '-' ;
    }
    return  max ;
}



