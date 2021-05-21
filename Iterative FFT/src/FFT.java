import java.util.Arrays;
import java.lang.*;

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

    private static int log2(int x) {
        //this code was taken from https://www.techiedelight.com/calculate-log-base-2-in-java/
        return (int) (Math.log(x) / Math.log(2));
    }

    private static int bitShuffle(int original, int n){
        int bits = log2(n);

        //convert to binary
        String binary = Integer.toBinaryString(original);

        //add leading zeros
        while(binary.length() != bits){
            binary = "0" + binary;
        }

        //reverse the string
        String reversed = "";

        for(int i = binary.length() - 1; i >= 0; i--)
        {
            reversed = reversed + binary.charAt(i);
        }

        //return integer version of the binary string created
        return Integer.parseInt(reversed,2);
    }

    private static Complex[] Iterative(int[] P){
        int n = P.length;
        int y = log2(n)+1;
        Complex[][] cache = new Complex[y][n];

        Complex[] Omega = getV(n*n,1);

        //initialize the cache
        for(int i = 0; i < n; i++){
            cache[0][i] = new Complex(bitShuffle(P[i], P.length),0);
        }

        for(int i = 1; i < y; i++){
            int size = (int) Math.pow(2,i);
//            size = size / 2;
            for(int j = 0; j < n  ; j += size){
                for(int k = 0; k < Math.floorDiv(size,2); k++){
                    Complex even = cache[i-1][j+k];
                    Complex odd = cache[i-1][j+k+Math.floorDiv(size,2)].times(Omega[(int) (Math.pow(2,i)*k)]);
                    cache[i][j+k] = even.plus(odd);
                    cache[i][j+k+Math.floorDiv(size,2)] = even.minus(odd);
                }
            }
        }
//        for(int i = 0; i < y; i++){
//            for(int j = 0; j < n; j++){
//                System.out.print(cache[i][j]+", ");
//            }
//            System.out.println();
//        }
        return cache[y-1];
    }

    private  static void compareFFTS(int[] poly, int n){
        boolean keepgoing = true;
        while(keepgoing){
            int[] numbers = new int[n];
            Complex[] test = new Complex[n];
            for(int i=0; i < n; i++){
                test[i] = new Complex(i,0);
                numbers[i] = i;
            }

            //calculate omega values
            Complex[] Omega = getV(n, 1);

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

            Complex[] IterativeAnswers = Iterative(numbers);

//        System.out.println("----------Forward FFT----------");
//        for(int i = 0; i < FFTAnswers.length; i++){
//            FFTAnswers[i] = FFTAnswers[i];
//            System.out.println(FFTAnswers[i]);
//        }

            endTime = System.nanoTime();
            long IterativeTime = (endTime - startTime)/1000000;


            if(FFTTime >= 300000|| IterativeTime >=300000){
                keepgoing = false;
            }

            System.out.println(n + "," + FFTTime + "," + IterativeTime);
            n= n*2;
        }
    }

    private static int[] makeArray(int PowerOf2){
        int size = (int) Math.pow(2,PowerOf2);
        int[] array = new int[size];
        for(int i=0; i < size; i++){
            array[i] = i;
        }
        return array;
    }

    public static void main(String[] args) {
        int power = 3;
        int[] numbers = makeArray(power);
        compareFFTS(numbers, (int) Math.pow(2,power));
//        Iterative(numbers);
    }
}