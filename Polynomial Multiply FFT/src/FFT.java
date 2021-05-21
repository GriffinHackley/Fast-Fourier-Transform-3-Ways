import java.util.Arrays;

public class FFT {

    private static Complex[][] split(Complex[] array){
        Complex[] even = new Complex[array.length/2];
        Complex[] odd = new Complex[array.length/2];

        for(int i=0; i < array.length; i++){
            if(i%2 == 0){
                even[i/2] = array[i];
            }
            else{
                odd[(i-1)/2] = array[i];
            }
        }
        Complex[][] split = new Complex[2][];
        split[0] = even;
        split[1] = odd;
        return split;
    }

    private static Complex[] getV(int n, int i){
        Complex[] Omega = new Complex[n];
        for(int j = 0; j < n; j++) {
            double theta = 2 * j * Math.PI / n;
            Omega[j] = new Complex(Math.cos(theta), Math.sin(i * theta));
        }
        return Omega;
    }

    private static Complex[] FFT(Complex[] P, Complex[] Omega, int n){
        if(n==1){
            return P;
        }

        //compute values for even and odd
        Complex[][] split = split(P);

        Complex[] even = FFT(split[0], Omega, n/2);
        Complex[] odd  = FFT(split[1], Omega, n/2);

        Omega = getV(n,1);

        //construct solution

        Complex[] answer = new Complex[n];

        for(int i = 0; i < n/2; i++){
            answer[i]       = even[i].plus(Omega[i].times(odd[i]));
            answer[i + n/2] = even[i].minus(Omega[i].times(odd[i]));
        }

        return answer;
    }

    private static Complex[] IFFT(Complex[] P, Complex[] Omega, int n){
        if(n==1){
            return P;
        }

        //compute values for even and odd
        Complex[][] split = split(P);

        Complex[] even = IFFT(split[0], Omega, n/2);
        Complex[] odd  = IFFT(split[1], Omega, n/2);


        //construct solution

        Omega = getV(n,-1);

        Complex[] answer = new Complex[n];

        for(int i = 0; i < n/2; i++){
            answer[i]       = even[i].plus(Omega[i].times(odd[i]));
            answer[i + n/2] = even[i].minus(Omega[i].times(odd[i]));
        }

        return answer;
    }

    private static double[] highSchool(double[] P, double[] Q){
        double[] answer = new double[P.length+Q.length];

        for(int i = 0; i < P.length; i++){
            for(int j = 0; j < Q.length; j++){
                answer[i+j] += P[i] * Q [j];
            }
        }
        return answer;
    }

    private static double[] threeSub(double[] P, double[] Q, int n){
        double[] answer = new double[(2*P.length)+1];

        if(n==1){
            answer[0] = P[0]*Q[0];
            return answer;
        }

        double[] PLow  = Arrays.copyOfRange(P,0,n/2);
        double[] PHigh = Arrays.copyOfRange(P,n/2, n);
        double[] QLow  = Arrays.copyOfRange(Q,0,n/2);
        double[] QHigh = Arrays.copyOfRange(Q,n/2, n);

        double[] tempP = new double[QHigh.length];
        double[] tempQ = new double[QHigh.length];

        for(int i = 0; i < QHigh.length; i++){
            tempP[i] = PLow[i] +  PHigh[i];
            tempQ[i] = QLow[i] +  QHigh[i];
        }

        double[] lowLow    = threeSub(PLow , QLow ,n/2);
        double[] highHigh  = threeSub(PHigh , QHigh,n/2);
        double[] sub       = threeSub(tempP , tempQ,n/2);

        double[] middle = new double[Q.length];

        for(int i = 0; i < Q.length; i++){
            middle[i] = sub[i]-lowLow[i]-highHigh[i];
        }


        for(int i = 0; i < n; i++){
            answer[i]       += lowLow[i];
            answer[i+(n/2)] += middle[i];
            answer[i+n]     += highHigh[i];
//            System.out.println("N: "+ n + " i: " + i);
        }
        return answer;
    }

    private static Complex[] FFTPolyMult(Complex[] P, Complex[] Q, int n){
        //run both arrays through FFT
        Complex[] Omega = getV(2*n, 1);
        Complex[] FFTPanswer = FFT(P,Omega,2*n);
        Complex[] FFTQanswer = FFT(P,Omega,2*n);

        //multiply answers together
        Complex[] FFTAnswer = new Complex[2*n];
        for(int i = 0; i < 2*n; i++){
            FFTAnswer[i] = FFTPanswer[i].times(FFTQanswer[i]);
        }

        //run answer through inverse FFT
        Complex[] OmegaInv = getV(n, -1);
        Complex[] IFFTAnswer = IFFT(FFTAnswer,OmegaInv, 2*n);

        //divide inverse FFT by 2*n
        for(int i=0; i < 2 * n; i++){
            IFFTAnswer[i] = IFFTAnswer[i].divides(new Complex(2*n,0));
        }
        return IFFTAnswer;
    }

    private static void runTest(){
        boolean keepgoing = true;
        int n = 32;
        while(keepgoing){
            Complex[] test = new Complex[n];
            for(int i=0; i < n; i++){
                test[i] = new Complex(i,0);
            }

            //calculate omega values
            Complex[] Omega = getV(n, 1);
            Complex[] OmegaInv = getV(n, -1);


            long startTime = System.nanoTime();

            Complex[] FFTAnswers = FFT(test, Omega, n);

//        System.out.println("----------Forward FFT----------");
//        for(int i = 0; i < FFTAnswers.length; i++){
//            FFTAnswers[i] = FFTAnswers[i];
//            System.out.println(FFTAnswers[i]);
//        }

            long endTime = System.nanoTime();
            long FFTTime = (endTime - startTime)/1000000;



            startTime = System.nanoTime();

            Complex[] answers = IFFT(FFTAnswers,OmegaInv, n);


            endTime = System.nanoTime();
            long IFFTTime = (endTime - startTime)/1000000;

            if(IFFTTime >= 300000 || FFTTime >= 300000){
                keepgoing = false;
            }

            System.out.println(n + "," + FFTTime + "," + IFFTTime);
            n= n*2;
        }
    }

    private static double AvgDiff (double[] correct, int[] incorrect){
        double total = 0;
        for(int i = 0; i<correct.length; i++ ){
            total += Math.abs(correct[i] - incorrect[i]);
        }
        return total/correct.length;
    }

    private static int testHighSchool(int minutes){
        int n = 1;
        long timeInSeconds = 0;
        long startTime = System.nanoTime();
        boolean keepGoing = true;
        while(keepGoing){
            double[] P = new double[n];
            double[] Q = new double[n];

            for(int i=0; i < n; i++){
                P[i] = i;
                Q[i] = i;
            }

            double[] answer = highSchool(P,Q);
            long endTime = System.nanoTime();
            long finalTime = (endTime - startTime)/1000000;
            timeInSeconds = finalTime/1000;
            System.out.println(n + ":" + timeInSeconds);

            if(timeInSeconds/60 >= minutes){
                keepGoing = false;
            }
            else{
                n*=2;
            }
        }
        System.out.println("t:" + timeInSeconds);
        return n;
    }

    private static int test3Sub(int minutes){
        int n = 1;
        long timeInSeconds = 0;
        long startTime = System.nanoTime();
        boolean keepGoing = true;
        while(keepGoing){
            double[] P = new double[n];
            double[] Q = new double[n];

            for(int i=0; i < n; i++){
                P[i] = i;
                Q[i] = i;
            }

            double[] answer = threeSub(P,Q,n);
            long endTime = System.nanoTime();
            long finalTime = (endTime - startTime)/1000000;
            timeInSeconds = finalTime/1000;
            System.out.println(n + ":" + timeInSeconds);

            if(timeInSeconds/60 >= minutes){
                keepGoing = false;
            }
            else{
                n*=2;
            }
        }
        System.out.println("t:" + timeInSeconds);
        return n;
    }

    private static int testFFT(int minutes){
        int n = 1;
        long timeInSeconds = 0;
        long startTime = System.nanoTime();
        boolean keepGoing = true;
        while(keepGoing){
            Complex[] FFTP = new Complex[2*n];
            Complex[] FFTQ = new Complex[2*n];

            //initialize polynomial arrays
            for(int i=0; i < n; i++){
                FFTP[i] = new Complex(i,0);
                FFTQ[i] = new Complex(i,0);
            }

            //pad the arrays that are going to be used in FFT with 0s
            for(int i = n; i < 2*n; i++){
                FFTP[i] = new Complex(0,0);
                FFTQ[i] = new Complex(0,0);
            }


            Complex[] answer= FFTPolyMult(FFTP,FFTQ,n);
            long endTime = System.nanoTime();
            long finalTime = (endTime - startTime)/1000000;
            timeInSeconds = finalTime/1000;

            if(timeInSeconds/60 >= minutes){
                keepGoing = false;
            }
            else{
                n*=2;
            }
            System.out.println(n + ":" + timeInSeconds);
        }
        System.out.println("t:" + timeInSeconds);
        return n;
    }

    private static void FFTVSHighSchool(){
        int n = 1;
        while(n < 10000){

            double[] P = new double[n];
            double[] Q = new double[n];

            Complex[] FFTP = new Complex[2*n];
            Complex[] FFTQ = new Complex[2*n];

            //initialize polynomial arrays
            for(int i=0; i < n; i++){
                FFTP[i] = new Complex(i,0);
                FFTQ[i] = new Complex(i,0);
                P[i] = i;
                Q[i] = i;
            }

            //pad the arrays that are going to be used in FFT with 0s
            for(int i = n; i < 2*n; i++){
                FFTP[i] = new Complex(0,0);
                FFTQ[i] = new Complex(0,0);
            }

            double[] highschoolAnswer = highSchool(P,Q);

            Complex[] FFTAnswer = FFTPolyMult(FFTP,FFTQ,n);

            //get only the real parts of the answer
            int[] FinalAnswer = new int[2*n];
            for(int i=0; i < 2*n; i++){
                FinalAnswer[i] = (int) FFTAnswer[i].re();
            }
            int[]casted = new int[highschoolAnswer.length];
            for(int i = 0; i < highschoolAnswer.length; i++){
                casted[i] = (int) highschoolAnswer[i];
            }

//            System.out.println("Highschool:" + Arrays.toString(casted));
//            System.out.println("FFT       :" + Arrays.toString(FinalAnswer));

            System.out.println(n + ": " + AvgDiff(highschoolAnswer,FinalAnswer));

            n*=2;
        }
    }

    private static void testAlg(int minutes){
        int threeSub = test3Sub(minutes);
        System.out.println("done with 3 sub");
        int HighSchool = testHighSchool(minutes);
        System.out.println("done with 3 highschool");
        int FFT = testFFT(minutes);
        System.out.println("done with FFT");

        System.out.println("Largest N gotten to in " + minutes + " minutes");
        System.out.println("3 sub      :" + threeSub/2);
        System.out.println("High School:" + HighSchool/2);
        System.out.println("FFT        :" + FFT/2);
    }

    public static void main(String[] args){
        FFTVSHighSchool();
//        testAlg(20);
    }
}