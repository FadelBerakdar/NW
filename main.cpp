/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  08/08/2016 16:34:43
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (),
 *   Organization:
 *
 * =====================================================================================
 */

#include "NW.h"

int main() {
    
    // TGCTAGTATAAACCTTATGGTATCTGCAGCAGAGGTTTCTTTAATCTCTCAATAGTAGATGCTTTGAAAC
    // TTATCTATAATTTGGTATTGTAATGACAGTTTGTGTTTGGTTTTTTCTTCAGTAT
    
    
    Alignment align("AACACT", "GCAA", 10, 1);
    align.global();
    
    Alignment a("TGCTAGTATAAACCTTATGGTATCTGCAGCAGAGGTTTCTTTAATCTCTCAATAGTAGATGCTTTGAAAC",
                "TTATCTATAATTTGGTATTGTAATGACAGTTTGTGTTTGGTTTTTTCTTCAGTAT", 10, 1);
    a.global();
    
}

