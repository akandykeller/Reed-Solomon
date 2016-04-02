//
//  EE/ME/CS 127 Final Project
//
//  Reed-Solomon Encoder and Decoder
//
//  Created by Thomas Andy Keller on 6/5/13.
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define REMAINDER 0
#define QUOTIENT 1

#define ELP 0
#define EEP 1

#define NONE -32

#define DECODING_ERROR -2

// Length of the code.
#define LENGTH 31

// alpha^5 as an integer in vector representation
#define A5 32

// alpha^3 + 1 as an integer in vector representation
#define A31 9

// Field elements are used in vector representation and stored as 32-bit integers.
typedef int FieldElement;

// Polynomials are stored as arrays of field elements. Each element of the array represents
// the coefficient of the corresponding x term. i.e.  if Polynomial[1] = alpha^2, this
// represents alpha^2 * x.
typedef FieldElement* Polynomial;

void printPoly(Polynomial poly);
FieldElement multiplyElements(FieldElement e1, FieldElement e2);


/*
 * Finite Field Element Operations
 */


//****************************************************************************************************
//
//  addElements(FieldElement e1, FieldElement e2);
//      Takes two field elements in vector representation and returns the sum in a binary field.
//
//****************************************************************************************************

FieldElement addElements(FieldElement e1, FieldElement e2) {
    
    // Return the bitwise xor'd vector representation of the sum if neither are zero
    if (e1 != NONE && e2 != NONE)
        return e1 ^ e2;
    
    // If one or both of the elements are zero, return non-zero element, or NONE
    else if (e1 == NONE && e2 != NONE)
        return e2;
    else if (e1 != NONE && e2 == NONE)
        return e1;
    else
        return NONE;
}


//****************************************************************************************************
//
//  toFieldElement(int aPower);
//     Takes an integer power of alpha and converts it to the corresponding 32-bit vector representing
//     a field element.
//
//****************************************************************************************************

FieldElement toFieldElement(int aPower) {
    FieldElement result = 1;
    
    // Catch where the alpha term is 0
    if (aPower == NONE)
        return 0;
    
    // Nonzero elements of the field form a multiplicitive cyclic group.
    while (aPower < 0) {
        aPower += LENGTH;
    }
    
    
    // To convert the power of alpha to a polynomial, let the polynomial equal 1, and then
    // continually multiply by alpha until alpha^5 is reached. Then remove it and add alpha^3 + 1.
    // Continue until no more powers of alpha remain.
    while (aPower > 0) {
        aPower--;
        result = result << 1;
        if ((result & A5) == A5) {
            result = addElements(result, A5);           // Addition = subtraction in binary field
            result = addElements(result, A31);
        }
    }
    
    return result;
}


//****************************************************************************************************
//
//  toAlphaPower(FieldElement e);
//     Takes 32-bit field element in vector form and returns the corresponding power of alpha.
//
//****************************************************************************************************

int toAlphaPower(FieldElement e) {
    int i;
    
    // Catch where the element is 0
    if (e == 0)
        return NONE;
    
    // Determine alpha power through use of already defined function
    for (i = 0; i < LENGTH; i++)
        if (toFieldElement(i) == e)
            return i;
    
    // Attempt negative values
    for (i = -LENGTH; i < 0; i++)
        if (toFieldElement(i) == e)
            return i;
    
    // Catch all
    return NONE;
}


//****************************************************************************************************
//
//  multiplyElements(FieldElement e1, FieldElement e2);
//      Takes two FieldElements as arguments, and returns their product in the F_32 field
//      with the primitive polynomial a^5 + a^3 + 1.
//
//****************************************************************************************************

FieldElement multiplyElements(FieldElement e1, FieldElement e2) {
    FieldElement product;
    int aPower1, aPower2, productPower;
    
    aPower1 = toAlphaPower(e1);
    aPower2 = toAlphaPower(e2);
    
    productPower = aPower1 + aPower2;
        
    // If either element is NONE, return NONE;
    if (e1 == NONE || e2 == NONE)
        return NONE;
    
    if (e1 == 0 || e2 == 0)
        return 0;
    
    
    if (e1 == 1)
        return e2;
    if (e2 == 1)
        return e1;
    

    product = toFieldElement(productPower);
    
    return product;
}


//****************************************************************************************************
//
//  divideElements(FieldElement e1, FieldElement e2);
//      Takes two FieldElements as arguments and returns e1 / e2 as a field element.
//
//****************************************************************************************************

FieldElement divideElements(FieldElement e1, FieldElement e2) {
    FieldElement result;
    
    int pow1 = toAlphaPower(e1);
    int pow2 = toAlphaPower(e2);
    
    result = toFieldElement(pow1 - pow2);
    
    return result;
}


//****************************************************************************************************
//
//  elementPower(FieldElement e, int power);
//      Takes a FieldElement and an integer power as arguments and returns the element raised to that
//      power as a FieldElement.
//
//****************************************************************************************************

FieldElement elementPower(FieldElement e, int power) {
    int i;
    FieldElement result = e;

    if (power == 0)
        return 1;
    
    // If element is NONE, return;
    if (e == NONE)
        return e;
    
    // Loop multiplication power times.
    for (i = 1; i < power; i++) {
        result = multiplyElements(result, e);
    }
    
    return result;
}


//****************************************************************************************************
//
//  randomElement();
//      returns a random non-zero FieldElement
//
//****************************************************************************************************

FieldElement randomElement() {
    FieldElement a = 0;
    
    // Return a non-zero element.
    while (a == 0)
        a = toFieldElement(rand() % 32);
    
    return a;
}


/*
 *  Polynomial Operations
 */


//****************************************************************************************************
//
//  initPolynomial();
//      Takes no arguments and returns a polynomial of length LENGTH with all elements equal to NONE.
//
//****************************************************************************************************

Polynomial initPolynomial() {
    int i;
    Polynomial poly = malloc(sizeof(FieldElement) * LENGTH);
    
    // Set all elements to NONE
    for (i = 0; i < LENGTH; i ++)
        poly[i] = NONE;
    
    return poly;
}

//****************************************************************************************************
//
//  polyDegree(Polynomial poly);
//      Takes a Polynomial and returns its integer degree. 
//
//****************************************************************************************************

int polyDegree(Polynomial poly) {
    int i = LENGTH - 1;
    
    // Find the first NONE item, and then search back to find the
    // first zero item.
    while (poly[i] == NONE || poly[i] == 0)
        i--;
    
    // decrement i before returning to give accurate degree.
    // NOTE: Returns -1 if all polynomail elements = NONE
    return i;
}

//****************************************************************************************************
//
//  isZero(Polynomial poly);
//      Takes a Polynomial and returns 1 if the polynomial is 0, else return 0.
//
//****************************************************************************************************

int isZero(Polynomial poly) {
    if(polyDegree(poly) <= 0 && poly[0] <= 0)
        return 1;
    
    return 0;
}


//****************************************************************************************************
//
//  addPolynomials(Polynomial poly1, Polynomial poly2);
//      Takes two Polynomials as arguments and returns their sum as a Polynomial.
//
//****************************************************************************************************

Polynomial addPolynomials(Polynomial poly1, Polynomial poly2) {
    int i;
    Polynomial polySum = initPolynomial();

    // Cycle through all elements and add 
    for (i = 0; i < LENGTH; i++)
        polySum[i] = addElements(poly1[i], poly2[i]);

    return polySum;
}


//****************************************************************************************************
//
//  multiplyPolynomials(Polynomial poly1, Polynomial poly2);
//      Takes two Polynomials as arguments and returns their product as a Polynomial.
//
//****************************************************************************************************

Polynomial multiplyPolynomials(Polynomial poly1, Polynomial poly2) {
    int i, j;
    
    int pow1 = polyDegree(poly1);
    int pow2 = polyDegree(poly2);
    
    Polynomial productFinal = initPolynomial();
    
    // Cycle through all elements of the first polynomial and distribute them through the
    // second polynomial.
    for (i = 0; i <= pow1; i++) {
        // Distribute element i of the first polynomial through the second polynomial and
        // store results in i + j = sum of corresponding x powers.
        
        // Set up empty temporary polynomial to hold partial product
        Polynomial productTemp = initPolynomial();
        
        for (j = 0; j <= pow2; j++) {
            productTemp[i + j] = multiplyElements(poly1[i], poly2[j]);
        }
        
        // Add the temporary product to the running sum
        productFinal = addPolynomials(productFinal, productTemp);
    }
    
    return productFinal;
}


//****************************************************************************************************
//
//  dividePolynomials(Polynomial poly1, Polynomial poly2, int returnType);
//      Takes two Polynomials and computes the quotient and remainder of poly1 / poly 2. Then, based
//      on the return setting it either returns the quotient, remainder or NULL.
//
//****************************************************************************************************

Polynomial dividePolynomials(Polynomial poly1, Polynomial poly2, int returnType) {
    int i;
    
    Polynomial remainder = poly1;
    Polynomial quotient = initPolynomial();
    
    // Cannot divide by zero.
    if (isZero(poly2)) {
        if (returnType == REMAINDER)
            return remainder;
        else if (returnType == QUOTIENT)
            return quotient;
        else
            return NULL;
    }
    
    /*
     * Following standard polynomial division algorithms, we set the remainder to be equal to 
     * the numerator, divide the leading element of this by the denominator, multiply the denominator
     * by this element, and then subtract this from remainder. This process is repeated until the 
     * remainder is zero or until the degree of the remainder drops below the degree of the denominator.
     */
    while (polyDegree(remainder) >= polyDegree(poly2) && !isZero(remainder)) {
        int degR = polyDegree(remainder);
        int degP2 = polyDegree(poly2);
        int degT = degR - degP2;
        
        // Element t holds the quantity which we multiply the denominator by at each step.
        // printf("remainder[degR] = a^%d, poly2[degP2] = a^%d \n", toAlphaPower(remainder[degR]), toAlphaPower(poly2[degP2]));
        FieldElement t = divideElements(remainder[degR], poly2[degP2]);
        // printf("remainder[degR]/poly2[degP2] = a^%d \n", toAlphaPower(t));
        
        // Store this element at the correct degree of x in the quotient. 
        quotient[degT] = addElements(quotient[degT], t);
        
        // printf("quotient : ");
        // printPoly(quotient);
        
        // Multiply the denominator by the element t.
        Polynomial k = initPolynomial();
        for (i = 0; i <= degP2; i++)
            // Each element is stored at degT + i to preserve powers of x.
            k[degT + i] = multiplyElements(poly2[i],t);
        
        // printf("k : ");
        // printPoly(k);
        
        // The remainder is then computed through the subtraction of these two polynomials.
        remainder = addPolynomials(remainder, k);
        // printf("remainder : ");
        // printPoly(remainder);
    }
    
    if (returnType == REMAINDER)
        return remainder;
    else if (returnType == QUOTIENT)
        return quotient;
    else
        return NULL;
    
}


//****************************************************************************************************
//
//  evaluatePolynomial(Polynomial poly, FieldElement e);
//      Takes a polynomial and evaluates it at x = e. It returns the FieldElement result.
//
//****************************************************************************************************

FieldElement evaluatePolynomial(Polynomial poly, FieldElement e) {
    int i;
    
    Polynomial result = initPolynomial();
    FieldElement eTemp = 0;
    FieldElement eResult = 0;
    
    // Plug in powers of alpha at x^i and then multiply by the polynomial coefficient
    for (i = 0; i < LENGTH; i++) {
        eTemp = elementPower(e, i);
        result[i] = multiplyElements(poly[i], eTemp);
        eResult = addElements(eResult, result[i]);
    }
                                 
    return eResult;
}


//****************************************************************************************************
//
//  findRoots(Polynomial poly);
//      Takes a polynomial and returns the roots.
//
//****************************************************************************************************

Polynomial findRoots(Polynomial poly) {
    int i;
    
    Polynomial roots = initPolynomial();
    
    // Plug in powers of alpha at x^i and then multiply by the polynomial coefficient
    for (i = 0; i < LENGTH; i++) {
        if (evaluatePolynomial(poly, toFieldElement(i)) == 0)
            roots[i] = toFieldElement(i);
    }
    
    return roots;
}


//****************************************************************************************************
//
//  formalDerivative(Polynomial poly);
//      Takes a polynomial and evaluates it at x = e. It returns the FieldElement result.
//
//****************************************************************************************************

Polynomial formalDerivative(Polynomial poly) {
    int i, j;
    int deg = polyDegree(poly);
    Polynomial result = initPolynomial();
    
    for (i = 0; i <= deg; i++) {
        for (j = 0; j <= i; j++)
            result[i] = addElements(result[i], poly[i + 1]);
    }
    
    return result;
}



/*
 * Code operations
 */

//****************************************************************************************************
//
//  encode(Polynomial message, int t);
//      Takes two Polynomials as arguments and returns their product as a Polynomial.
//
//****************************************************************************************************

Polynomial encode(Polynomial msg, int t) {
    int i;
    
    Polynomial codeword = initPolynomial();
    Polynomial mx = initPolynomial();
    Polynomial rmx = initPolynomial();
    Polynomial g = initPolynomial();
    Polynomial gTemp = initPolynomial();
    
    int xPow = 2*t;
    
    // multiply the message vector by x^(n-k)
    for (i = 0; i < LENGTH; i++) {
        mx[i + xPow] = msg[i];
    }
    
    // Set up g(x) with first binomial.
    g[0] = toFieldElement(1);
    g[1] = toFieldElement(0);
    
    // Multiply subsequent binomials until final generator is achieved. 
    for (i = 2; i <= 2*t; i++) {
        gTemp[0] = toFieldElement(i);
        gTemp[1] = toFieldElement(0);
        
        g = multiplyPolynomials(g, gTemp);
    }
    
    // get the remainder after division.
    rmx = dividePolynomials(mx, g, REMAINDER);
    
    codeword = addPolynomials(mx, rmx);
    
    return codeword;
}


//****************************************************************************************************
//
//  syndrome(Polynomial recieved);
//      Takes a recieved polynomial as an argument and returns the syndrome polynomial.
//
//****************************************************************************************************

Polynomial syndrome(Polynomial recieved, int t) {
    int i;

    Polynomial synd = initPolynomial();
    FieldElement e;
    
    for (i = 0; i < 2*t; i++) {
        e = toFieldElement(i + 1);
        synd[i] = evaluatePolynomial(recieved, e);
    }
    
    return synd;
}



//****************************************************************************************************
//
//  euclideanAlgorithm(Polynomial recieved, int t, int returnType);
//      Takes a syndrome, the error correcting capability of the code, and computes the error
//      locator polynomial as well as omega. Then, based on the returnType argument, it returns the
//      desired polynomial.
//
//****************************************************************************************************

Polynomial euclideanAlgorithm(Polynomial recieved, int t_int, int returnType) {
    int xPow = 2*t_int;
    
    // Set up a ton of polynomials to hold all the incremental values used in the recursive algorithm.
    Polynomial quotient = initPolynomial();
    Polynomial remainder = initPolynomial();
    Polynomial s = initPolynomial();
    Polynomial t = initPolynomial();
    Polynomial s_0 = initPolynomial();
    Polynomial t_0 = initPolynomial();
    Polynomial s_N1 = initPolynomial();
    Polynomial t_N1 = initPolynomial();
    Polynomial b = initPolynomial();
    Polynomial a = initPolynomial();

    // Initialize the syndrome and a
    b = syndrome(recieved, t_int);
    
    remainder = b;
    
    a[xPow] = 1;
    
    // Set up the initial values for the recusrive extended Euclidean algorithm
    s_0[0] = 0;
    s_N1[0] = 1;
    t_0[0] = 1;
    t_N1[0] = 0;
    
    t = t_0;
    
    while (polyDegree(remainder) >= t_int) {
        
        quotient = dividePolynomials(a, b, QUOTIENT);
        remainder = addPolynomials(a, multiplyPolynomials(b, quotient));     //dividePolynomials(a, b, REMAINDER);
        s = addPolynomials(s_N1, multiplyPolynomials(quotient, s_0));
        t = addPolynomials(t_N1, multiplyPolynomials(quotient, t_0));
        
        // increment the s and t values to allow for recursive functionality.
        s_N1 = s_0;
        s_0 = s;
        t_N1 = t_0;
        t_0 = t;
        
        // increment remainder
        a = b;
        b = remainder;
    }
    
    if (returnType == ELP)
        return t;
    else if (returnType == EEP)
        return remainder;
    else
        return NULL;
}


//****************************************************************************************************
//
//  sendThroughChannel(Polynomial codeword, int p);
//      Simulates sending a given codeword through our channel model where symbol errors occur
//      independantly with probability p/100. Returns the transmitted codeword + random errors.
//
//****************************************************************************************************

Polynomial sendThroughChannel(Polynomial recieved, int p) {
    int i;
    FieldElement a;
    Polynomial transmitted = initPolynomial();
        
    // Loop through each element and add random error if probability is right.
    for (i = 0; i < LENGTH; i++) {
        if ((rand() % 100) < p)
            a = randomElement();
        else
            a = 0;
        
        transmitted[i] = addElements(recieved[i], a);
    }
    
    return transmitted;
}

//****************************************************************************************************
//
//  getErrorIndcies(Polynomial roots, int* numErrors, int* indicies, int* rootIndicies);
//      Takes the roots polynomial and a pointer to an integer counter and returns an array filled
//      with the indicies to the various errors. Additionally it fills the integer counter with the
//      number of errors.
//
//****************************************************************************************************

void getErrorIndicies(Polynomial roots, int* numErrors, int* indicies, int* rootIndicies) {
    int i;
    int j = 0;
    
    for (i = 0; i < LENGTH; i++) {
        if (roots[i] != 0 && roots[i] != NONE) {
            indicies[j] = toAlphaPower(toFieldElement(-i));
            rootIndicies[j] = i;
            j++;
        }
    }
    
    *numErrors = j;
}


//****************************************************************************************************
//
//  errorValues(Polynomial elp, Polynomial eep, Polynomial roots);
//      Takes the error locator polynomial and the error evaluate polynomial, and evaluates them using
//      forney's algorithm to compute the error values.
//
//****************************************************************************************************

Polynomial errorValues(Polynomial elp, Polynomial eep, Polynomial roots) {
    int i;
    
    Polynomial elp_prime = initPolynomial();
    Polynomial errors = initPolynomial();
    
    int errorIndicies[LENGTH] = {0};
    int rootIndicies[LENGTH] = {0};
    int numErrors;
    
    getErrorIndicies(roots, &numErrors, errorIndicies, rootIndicies);
    
    FieldElement elp_prime_e;
    FieldElement eep_e;

    elp_prime = formalDerivative(elp);
    
    for (i = 0; i < numErrors; i++) {
        elp_prime_e = evaluatePolynomial(elp_prime, roots[rootIndicies[i]]);
        
        if (elp_prime_e == 0) {
            //printf("Decoding Error \n");
            errors[errorIndicies[i]] = DECODING_ERROR;
        }
        else {
            eep_e = evaluatePolynomial(eep, roots[rootIndicies[i]]);
            errors[errorIndicies[i]] = divideElements(eep_e, elp_prime_e);
        }
    }
    
    return errors;
}


//****************************************************************************************************
//
//  decode(Polynomial recieved, int t);
//      Takes a recieved message and the error correcting capability, t, of the code, and returns the
//      decoded vector or NULL if decoding failed. This function uses the extended Euclidean
//      algorithm to complete decoding.
//
//****************************************************************************************************

Polynomial decode(Polynomial recieved, int t) {
    int i;
    int numErrors = 0;
    int errorIndicies[LENGTH] = {0};
    int rootIndicies[LENGTH] = {0};
    
    Polynomial decoded = recieved;
    Polynomial elp = initPolynomial();
    Polynomial eep = initPolynomial();
    Polynomial roots = initPolynomial();
    Polynomial errors = initPolynomial();
    
    
    
    elp = euclideanAlgorithm(recieved, t, ELP);
    
    roots = findRoots(elp);

    eep = euclideanAlgorithm(recieved, t, EEP);
    
    errors = errorValues(elp, eep, roots);
    
    getErrorIndicies(roots, &numErrors, errorIndicies, rootIndicies);
    
    for (i = 0; i < numErrors; i++) {
        if (errors[errorIndicies[i]] == DECODING_ERROR) {
            return NULL;
        }
        decoded[errorIndicies[i]] = addElements(decoded[errorIndicies[i]], errors[errorIndicies[i]]);
    }
    
    return decoded;
}

//****************************************************************************************************
//
//  getDecodedMessage(Polynomial recieved, int t);
//      Takes the decoded codeword and returns the last k elements which are the encoded message
//      given that the encoder is systematic.
//
//****************************************************************************************************

Polynomial getDecodedMessage(Polynomial decoded, int t) {
    int i;
    
    Polynomial decoded_msg = initPolynomial();
    
    for (i = 0; i < LENGTH - 2 * t; i++) {
        decoded_msg[i] = decoded[2 * t + i];
    }
    
    return decoded_msg;
}

//****************************************************************************************************
//
//  isEqual(Polynomial poly1, Polynomial poly2);
//      Takes two polynomials and returns 1 if they are equivalent, else 0.
//
//****************************************************************************************************

int isEqual(Polynomial poly1, Polynomial poly2) {
    int i;
   
    for (i = 0; i < LENGTH; i++) {
        if (poly1[i] != poly2[i]) {
            if ((poly1[i] == NONE && poly2 == 0) || (poly2[i] == NONE && poly1 == 0))
                continue;
            else
                return 0;
        }
    }
    
    return 1;
}



/*
 *  Testing and Simulation
 */

//****************************************************************************************************
//
//  toString(FieldElement e);
//      Takes a field element as an argument and returns its string representation.
//
//****************************************************************************************************

char* toString(FieldElement e) {
    
    char* str = malloc(sizeof(char) * 50);
    int first = 1;
    
    if ((e & 16) == 16) {
        strcat(str, "a^4");
        first = 0;
    }
    
    if ((e & 8) == 8) {
        if (!first)
            strcat(str, " + ");
        strcat(str, "a^3");
        first = 0;
    }
    
    if ((e & 4) == 4) {
        if (first != 1)
            strcat(str, " + ");
        strcat(str, "a^2");
        first = 0;
    }
    
    if ((e & 2) == 2) {
        if (!first)
            strcat(str, " + ");
        strcat(str, "a^1");
        first = 0;
    }
    
    if ((e & 1) == 1) {
        if (!first)
            strcat(str, " + ");
        strcat(str, "1");
    }
    
    if (e == 0)
        return "0";

    return str;
}


//****************************************************************************************************
//
//  printPoly(Polynomial poly1, int deg1);
//      Takes a polynomial and its power as arguments and prints a string representation of the
//      polynomial. Mainly for debugging.
//
//****************************************************************************************************

void printPoly(Polynomial poly) {
    int i;
    
    // Cycle through all elements of the polynomial and convert them to strings.
    if (poly[0] != NONE && poly[0] != 0)
        printf("(%s)  ", toString(poly[0]));
    else
        printf("0  ");
    for (i = 1; i <= polyDegree(poly); i++) {
        if (poly[i] != NONE && poly[i] != 0)
            printf("+  (%s)x^%d  ", toString(poly[i]), i);
    }
    printf("\n");
    
}

//****************************************************************************************************
//
//  simulation(int t, int p);
//      Runs a simulation to gather data on the decoding error rate as a function of p for a given t.
//      Returns 0 if decoding successful, returns 1 if there is a decoding error.
//
//****************************************************************************************************

int simulation(int t, int p) {
    
    Polynomial msg = initPolynomial();
    Polynomial codeword = initPolynomial();
    Polynomial transmitted = initPolynomial();
    Polynomial decoded = initPolynomial();
    Polynomial decoded_msg = initPolynomial();
    
    int i;
    
    // Generate random binary message vector of length k.
    for (i = 0; i < LENGTH - 2 * t; i++)
        msg[i] = rand() % 2;
    
    // Encode the message vector
    codeword = encode(msg,t);
    
    // Send it through the simulated channel model with a probability of error p/100.
    transmitted = sendThroughChannel(codeword, p);
    
    // Decode the transmitted code.
    decoded = decode(transmitted, t);
    
    // If NULL it means there was a decoding failure. 
    if (decoded == NULL)
        return 1;
        
    // Extract the encoded message.
    decoded_msg = getDecodedMessage(decoded, t);
    
    // Check if message was decoded incorrectly.
    if(!isEqual(msg, decoded_msg))
        return 1;
    
    // Else the decoding was correct
    else
        return 0;
    
}

//****************************************************************************************************
//
//  getErrorRate(int t, int p);
//      Given an error correcting value t, and the probability of error through a simulated channel p,
//      return the empirical error rate as a double. The simulation is run long enough for 100 errors
//      to happen.
//
//****************************************************************************************************

double getErrorRate(int t, int p) {
    
    int runCount = 0;
    int  error;
    double errorCount = 0.0;
    double errorRate;
    
    
    while (errorCount < 100.0 && runCount < 100000) {
        //printf("Run %d: ", runCount + 1);
        
        error = simulation(t, p);
        
        if (error) {
            //printf("Decoding Error\n");
            errorCount++;
        }
        else {
            //printf("Decoding Success\n");
        }
        
        runCount++;
    }
    
    errorRate = errorCount / runCount;
    
    printf("Error Rate = %f \n", errorRate);
    
    return errorRate;
}

//****************************************************************************************************
//
//  factorial(int n);
//      Takes an integer and returns its factorial.
//
//****************************************************************************************************

long double factorial(int n) {
    if (n == 1 || n == 0)
        return 1.0;
    else
        return (long double) n * factorial(n - 1);
}

//****************************************************************************************************
//
//  nChooseK(int n, int k);
//      Takes two integers and returns their equivalent binomial coefficient
//
//****************************************************************************************************

double nChooseK(int n, int k) {
    long double num = (double) factorial(n);
    long double denom = (double) (factorial(k) * factorial(n-k));
    
    
    return num/denom;
}

//****************************************************************************************************
//
//  getPredictedErrorRate(int t, int p);
//      Given an error correcting value t, and the probability of error through a simulated channel p,
//      return the theoretical error rate as a double.
//
//****************************************************************************************************

double getPredictedErrorRate(int t, int p) {
    
    double predictedErrorRate = 0.0;
    double P = (double) p/100;
    double a,b,c;
    double i;
    
    
    for (i = (double) t + 1.0; i < (double) LENGTH; i += 1.0) {
        a = nChooseK(LENGTH, i);
        b = pow(P, i);
        c = pow((1 - P), LENGTH - i);
       predictedErrorRate += a * b * c;
    }
    
    return predictedErrorRate;
}

//****************************************************************************************************
//
//  main(int argc, const char * argv[]);
//      Entry point of the program for testing and simulation.
//
//****************************************************************************************************

int main(int argc, const char * argv[]) {
    srand((int) time(NULL));
    int p, t;
    int i, j;

    double errorRates[3][30] = {0};
    double predictedErrorRates[3][30] = {0};

    // Calculate error rates for t = 1, 2, 3 from P = 0.01 to P = 0.30
    for (t = 1; t <= 3; t++)
        for (p = 0; p < 30; p += 1) {
            printf("Calculating error rate for t = %d, p = %d/100 \n", t, p + 1);
            errorRates[t - 1][p] = getErrorRate(t, p + 1);
        }
    
    
    // Print the computed error rates
    printf("Error Rates:\n");
    printf("p , rate\n");
    
    for (i = 0; i < 3; i++) {
        printf("t = %d\n", i + 1);
        for (j = 0; j < 30; j++) {
            printf("%d, %f\n", j + 1, errorRates[i][j]);
        }
    }
    
    // Compute the predicted error rates for the same values
    for (t = 1; t <= 3; t++)
        for (p = 0; p < 30; p += 1) {
            predictedErrorRates[t - 1][p] = getPredictedErrorRate(t, p + 1);
    }
    
    // Print them out.
    printf("Predicted Error Rates:\n");
    printf("p , rate\n");
    
    for (i = 0; i < 3; i++) {
        printf("t = %d\n", i + 1);
        for (j = 0; j < 30; j++) {
            printf("%d, %f\n", j + 1, predictedErrorRates[i][j]);
        }
    }
    
    
    return 0;
}


