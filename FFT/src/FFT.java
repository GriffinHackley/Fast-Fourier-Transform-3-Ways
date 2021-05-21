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


        //construct solution

        Complex[] answer = new Complex[n];

        for(int i = 0; i < n/2; i++){
            Omega = getV(n,1);
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

        Complex[] answer = new Complex[n];

        for(int i = 0; i < n/2; i++){
            Omega = getV(n,-1);
            answer[i]       = even[i].plus(Omega[i].times(odd[i]));
            answer[i + n/2] = even[i].minus(Omega[i].times(odd[i]));
        }
        return answer;
    }

    public static void main(String[] args) {
        boolean keepgoing = true;
        int n = 128;
        while(keepgoing){
            Complex[] test = new Complex[n];
            for(int i=0; i < n; i++){
                test[i] = new Complex(i,0);
            }

            //calculate omega values
            Complex[] Omega = getV(n, 1);
            Complex[] OmegaInv = getV(n, 1);


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

//        System.out.println();
//        System.out.println("----------Inverse FFT----------");
//        for(int i = 0; i < answers.length; i++){
//            answers[i] = answers[i].divides(new Complex(n,0));
//            System.out.println(answers[i]);
//        }

            endTime = System.nanoTime();
            long IFFTTime = (endTime - startTime)/1000000;

            if(IFFTTime >= 300000 || FFTTime >= 300000){
                keepgoing = false;
            }

            System.out.println(n + "," + FFTTime + "," + IFFTTime);
            n= n*2;
        }
    }
}