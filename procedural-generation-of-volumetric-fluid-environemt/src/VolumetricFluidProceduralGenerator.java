import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;

public class VolumetricFluidProceduralGenerator {
    public static void main(String[] args) throws IOException {
        System.out.println("Starting program.");

        // Width, height and depts of the volume
        int size = 128;
        // Height of the ground and height of the waves
        int groundPosition = 10;
        int waveHeight = 5;
        // This many volumes will be appended consecutively in the output file
        int simulationIterations = 10;

        VolumeGenerator vGenerator = new VolumeGenerator(size, groundPosition, waveHeight, simulationIterations);
        vGenerator.generateVolume();
    }
}

class VolumeGenerator {
    public int size;
    public float[][][] voxArray;
    public float[][][] voxOutputArray;
    public float[][][] voxBuffer;
    public float[][][] uVelocityArray;
    public float[][][] vVelocityArray;
    public float[][][] wVelocityArray;

    public float airPressure = 10f;
    public float waterPressure = 125f;
    public float groundPressure = 255f;
    public int groundPosition;
    public int waveHeight;
    public int simulationIterations;

    public VolumeGenerator(int size, int groundPosition, int waveHeight, int simulationIterations) {
        this.size = size;
        this.groundPosition = groundPosition;
        this.waveHeight = waveHeight;
        this.simulationIterations = simulationIterations;
        voxArray = new float[size][size][size];
        voxOutputArray = new float[size][size][size];
        voxBuffer = new float[size][size][size];
        uVelocityArray = new float[size][size][size];
        vVelocityArray = new float[size][size][size];
        wVelocityArray = new float[size][size][size];
    }

    // Main method of class - executes all other methods in the intended order
    public void generateVolume() throws IOException {

        generateVelocities();
        deformWater();

        boolean simulating = true;
        int count = 0;

        // Algorithm roughly modeled from Stams methods
        while(simulating){

            updateBuffer();
            diffuse(size, 0, voxArray, voxBuffer, 1f, 0.1f);

            updateBuffer();
            advect(size, voxArray, voxBuffer, uVelocityArray, vVelocityArray, wVelocityArray, 0.1f);

            updateOutputArray();
            overwriteFloor(groundPosition);
            overwriteAir(waveHeight);

            appendToOutputFile(simulationIterations, count);

            count++;
            if(count == simulationIterations)
                simulating = false;

        }
        createOutputFile(99);
    }

    // Copy voxArray into voxBuffer
    private void updateBuffer(){
        for(int x = 0; x < size; x++){
            for(int y = 0; y < size; y++){
                for(int z = 0; z < size; z++){
                    voxBuffer[x][y][z] = voxArray[x][y][z];
                }
            }
        }
    }

    // Copy voxArray into voxOutputArray
    private void updateOutputArray(){
        for(int x = 0; x < size; x++){
            for(int y = 0; y < size; y++){
                for(int z = 0; z < size; z++){
                    voxOutputArray[x][y][z] = voxArray[x][y][z];

                }
            }
        }
    }

    // J Stam methods:

    private void diffuse(int N, int b, float[][][] newArray, float[][][] oldArray, float dt, float diff){
        int nMax = N - 1;
        float a = dt*diff*nMax*nMax;
        for(int k = 0; k < 20; k++){
            for(int x = 1; x < nMax; x++){
                for(int y = 1; y < nMax; y++){
                    for(int z = 1; z < nMax; z++){
                        newArray[x][y][z] = (oldArray[x][y][z] + a *
                                (newArray[x - 1][y][z] + newArray[x + 1][y][z] +
                                newArray[x][y - 1][z] + newArray[x][y + 1][z] +
                                newArray[x][y][z - 1] + newArray[x][y][z + 1]))
                                / (1 + 6 * a);

                    }
                }
            }
            set_bnd(N);
        }
    }

    private void advect(int N, float[][][] newArray, float[][][] oldArray, float[][][] u, float[][][] v, float[][][] w, float dt){
        float x, y, z, s0, t0, s1, t1, r0, r1;
        int i0, j0, i1, j1, k0, k1;
        float dt0 = dt * N;

        for(int i = 0; i < N - 1; i++){
            for(int j = 0; j < N - 1; j++){
                for(int k = 0; k < N - 1; k++){
                    x = i - dt0 * u[i][j][k];
                    y = j - dt0 * v[i][j][k];
                    z = k - dt0 * w[i][j][k];

                    if(x < 0.5)
                        x = 0.5f;
                    if(x > N + 0.5){
                        x = N + 0.5f;
                    }
                    i0 = (int) x;
                    if(i0 >= 127) i0 = 126;
                    i1 = i0 + 1;

                    if(y < 0.5)
                        y = 0.5f;
                    if(y > N + 0.5){
                        y = N + 0.5f;
                    }
                    j0 = (int) y;
                    if(j0 >= 127) j0 = 126;
                    j1 = j0 + 1;

                    if(z < 0.5)
                        z = 0.5f;
                    if(z > N + 0.5){
                        z = N + 0.5f;
                    }
                    k0 = (int) z;
                    if(k0 >= 127) k0 = 126;
                    k1 = k0 + 1;

                    s1 = x - i0;
                    s0 = 1 - s1;
                    t1 = y - j0;
                    t0 = 1 - t1;
                    r1 = z - k0;
                    r0 = 1 - r1;

                    newArray[i][j][k] = s0 * t0 * r0 * oldArray[i0][j0][k0] +
                                        s0 * t1 * r0 * oldArray[i0][j1][k0] +
                                        s1 * t0 * r0 * oldArray[i1][j0][k0] +
                                        s1 * t1 * r0 * oldArray[i1][j1][k0] +
                                        s0 * t0 * r1 * oldArray[i0][j0][k1] +
                                        s0 * t1 * r1 * oldArray[i0][j1][k1] +
                                        s1 * t0 * r1 * oldArray[i1][j0][k1] +
                                        s1 * t1 * r1 * oldArray[i1][j1][k1];
                }
            }
        }
        set_bnd(N);
    }

    private void set_bnd(int N){
        int nMin = 0;
        int nMax = N - 1;

        // Popravi robove
        for(int i = 1; i <= nMax; i++){
            for(int j = 1; j <= nMax; j++){
                voxArray[i][j][0] = voxArray[i][j][1];
                voxArray[i][j][nMax] = voxArray[i][j][nMax - 1];

                voxArray[i][0][j] = voxArray[i][1][j];
                voxArray[i][nMax][j] = voxArray[i][nMax - 1][j];

                voxArray[0][i][j] = voxArray[1][i][j];
                voxArray[nMax][i][j] = voxArray[nMax - 1][i][j];
            }
        }

        // Popravi ogljišča
        voxArray[0][0][0] = 1/3f * (voxArray[1][0][0] + voxArray[0][1][0] + voxArray[0][0][1]);

        voxArray[nMax][0][0] = 1/3f * (voxArray[nMax - 1][0][0] + voxArray[nMax][1][0] + voxArray[nMax][0][1]);
        voxArray[0][nMax][0] = 1/3f * (voxArray[1][nMax][0] + voxArray[0][nMax - 1][0] + voxArray[0][nMax][1]);
        voxArray[0][0][nMax] = 1/3f * (voxArray[1][0][nMax] + voxArray[0][1][nMax] + voxArray[0][0][nMax - 1]);

        voxArray[nMax][nMax][0] = 1/3f * (voxArray[nMax - 1][nMax][0] + voxArray[nMax][nMax - 1][0] + voxArray[nMax][nMax][1]);
        voxArray[0][nMax][nMax] = 1/3f * (voxArray[1][nMax][nMax] + voxArray[0][nMax - 1][nMax] + voxArray[0][nMax][nMax - 1]);
        voxArray[nMax][0][nMax] = 1/3f * (voxArray[nMax - 1][0][nMax] + voxArray[nMax][1][nMax] + voxArray[nMax][0][nMax - 1]);
        voxArray[nMax][0][nMax] = 1/3f * (voxArray[nMax - 1][nMax][nMax] + voxArray[nMax][nMax - 1][nMax] + voxArray[nMax][nMax][nMax - 1]);
    }


    // Deform water with perlin noise
    private void deformWater(){
        PerlinNoiseGenerator fn = new PerlinNoiseGenerator(156437);
        // fn.SetNoiseType(FastNoiseLite.NoiseType.Perlin);
        for(int x = 0; x < size; x++){
            for(int y = 0; y < size; y++){
                for(int z = 0; z < size; z++){
                    float perlinNoiseFactor = fn.noise3(x * 0.01f, y * 0.01f, z * 0.01f);

                    voxArray[x][y][z] = waterPressure + perlinNoiseFactor * 75;

                }
            }
        }
    }

    // Generate velocitiet with perlin noise
    private void generateVelocities(){
        PerlinNoiseGenerator fnu = new PerlinNoiseGenerator(567547);
        PerlinNoiseGenerator fnv = new PerlinNoiseGenerator(345634);
        PerlinNoiseGenerator fnw = new PerlinNoiseGenerator(456345);

        float[][][] velocityArray = new float[size][size][size];
        float vMin = 1000;
        float vMax = -1000;
        float uMin = 1000;
        float uMax = -1000;
        float wMin = 1000;
        float wMax = -1000;
        for(int x = 0; x < size; x++){
            for(int y = 0; y < size; y++){
                for(int z = 0; z < size; z++){
                    float perlinNoiseFactorU = fnu.noise3(x * 0.01f, y * 0.01f, z * 0.01f);
                    float perlinNoiseFactorV = fnv.noise3(x * 0.01f, y * 0.01f, z * 0.01f);
                    float perlinNoiseFactorW = fnw.noise3(x * 0.01f, y * 0.01f, z * 0.01f);

                    uVelocityArray[x][y][z] = perlinNoiseFactorU;
                    vVelocityArray[x][y][z] = perlinNoiseFactorV;
                    wVelocityArray[x][y][z] = perlinNoiseFactorW;

                    if(uVelocityArray[x][y][z] > uMax) uMax = uVelocityArray[x][y][z];
                    if(uVelocityArray[x][y][z] < uMin) uMin = uVelocityArray[x][y][z];
                    if(vVelocityArray[x][y][z] > vMax) vMax = uVelocityArray[x][y][z];
                    if(vVelocityArray[x][y][z] < vMin) vMin = uVelocityArray[x][y][z];
                    if(wVelocityArray[x][y][z] > wMax) wMax = wVelocityArray[x][y][z];
                    if(wVelocityArray[x][y][z] < wMin) wMin = wVelocityArray[x][y][z];
                }
            }
        }
        float min = 100;
        float max = 0;
        System.out.println("vMax vMin diff: " + wMax + " " +  wMin);
        for(int x = 0; x < size; x++){
            for(int y = 0; y < size; y++){
                for(int z = 0; z < size; z++){
                    float uOffset = -1 * uMin;
                    if(uMin < 0)
                        uOffset = -1 * uMin;
                    if(uMin > 1)
                        uOffset = uMin;
                    uVelocityArray[x][y][z] += uOffset;
                    uVelocityArray[x][y][z] *= (1/ Math.abs(uMax - uMin));
                    //uVelocityArray[x][y][z] -= 0.5f;
                    float vOffset = -1 * vMin;
                    if(vMin < 0)
                        vOffset = -1 * vMin;
                    if(vMin > 1)
                        vOffset = vMin;
                    vVelocityArray[x][y][z] += vOffset;
                    vVelocityArray[x][y][z] *= (1/ Math.abs(vMax - vMin));
                    //vVelocityArray[x][y][z] -= 0.5f;
                    float wOffset = -1 * wMin;
                    if(wMin < 0)
                        wOffset = -1 * wMin;
                    if(wMin > 1)
                        wOffset = wMin;
                    wVelocityArray[x][y][z] += wOffset;
                    wVelocityArray[x][y][z] *= (1/ Math.abs(wMax - wMin));
                    //wVelocityArray[x][y][z] -= 0.5f;
                    if(wVelocityArray[x][y][z] > max) max = wVelocityArray[x][y][z];
                    if(wVelocityArray[x][y][z] < min) min = wVelocityArray[x][y][z];

                }
            }
        }
        System.out.println("Max min diff: " + max + " " +  min);

        /*
        for(int x = 0; x < size; x++){
            for(int y = 0; y < size; y++){
                for(int z = 0; z < size - 1; z++){

                    uVelocityArray[x][y][z] = velocityArray[x][y][z];// - velocityArray[x][y][z + 1];
                    vVelocityArray[x][y][z] = velocityArray[x][z][y];// - velocityArray[x][z + 1][y];
                    wVelocityArray[x][y][z] = velocityArray[z][x][y];// - velocityArray[z + 1][x][y];

                    if(z == size - 2){
                        uVelocityArray[x][y][z + 1] = 0f;
                        vVelocityArray[x][y][z + 1] = 0f;
                        wVelocityArray[x][y][z + 1] = 0f;
                    }
                }
            }
        }*/
        // System.out.println(uVelocityArray[60][60][60]);

    }

    // Calculate gradients with code from script
    private float[] computeGradientsScript(int x, int y, int z){

        if(x == 0 || x == size - 1 || y == 0 || y == size - 1 || z == 0 || z == size - 1)
            return new float[]{0f,0f,0f};

        float[] kernelX = new float[]{
                1, 0, -1,
                2, 0, -2,
                1, 0, -1,

                2, 0, -2,
                4, 0, -4,
                2, 0, -2,

                1, 0, -1,
                2, 0, -2,
                1, 0, -1
        };
        float[] kernelY = new float[]{
                1,  2,  1,
                0,  0,  0,
                -1, -2, -1,

                2,  4,  2,
                0,  0,  0,
                -2, -4, -2,

                1,  2,  1,
                0,  0,  0,
                -1, -2, -1
        };
        float[] kernelZ = new float[]{
                1,  2,  1,
                2,  4,  2,
                1,  2,  1,

                0,  0,  0,
                0,  0,  0,
                0,  0,  0,

                -1, -2, -1,
                -2, -4, -2,
                -1, -2, -1
        };
        float[] values = new float[]{
                0f, 0f, 0f,
                0f, 0f, 0f,
                0f, 0f, 0f,

                0f, 0f, 0f,
                0f, 0f, 0f,
                0f, 0f, 0f,

                0f, 0f, 0f,
                0f, 0f, 0f,
                0f, 0f, 0f
        };

        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                for(int k = 0; k < 3; k++){
                    int index = i + j * 3 + k * 3 * 3;
                    values[index] = voxOutputArray[x + i - 1][y + j - 1][z + k - 1] / 255;
                }
            }
        }

        float dx = 0f;
        float dy = 0f;
        float dz = 0f;
        int len = kernelX.length;
        for (int i = 0; i < len; i++) {
            dx += kernelX[i] * values[i];
            dy += kernelY[i] * values[i];
            dz += kernelZ[i] * values[i];
        }
        dx /= len;
        dy /= len;
        dz /= len;
        dx *= 100;
        dy *= 100;
        dz *= 100;
        return new float[]{dx, dy, dz};
    }

    // Calculate gradients with code from SO - doesn't work
    private float[] computeGradients(int x, int y, int z){
        float[] res = new float[3];
        if(x == 0 || x == size - 1 || y == 0 || y == size - 1 || z == 0 || z == size - 1)
            return new float[]{0f,0f,0f};
        // x
        int xi = (int)(x + 0.5f);
        float xf = x + 0.5f - xi;
        float xd0 = voxOutputArray[xi - 1][y][z];
        float xd1 = voxOutputArray[xi][y][z];
        float xd2 = voxOutputArray[xi + 1][y][z];
        res[0] = (xd1 - xd0) * (1.0f - xf) + (xd2 - xd1) * xf; // lerp
        // y
        int yi = (int)(y + 0.5f);
        float yf = y + 0.5f - yi;
        float yd0 = voxOutputArray[x][yi - 1][z];
        float yd1 = voxOutputArray[x][yi][z];
        float yd2 = voxOutputArray[x][yi + 1][z];
        res[1] = (yd1 - yd0) * (1.0f - yf) + (yd2 - yd1) * yf; // lerp
        // z
        int zi = (int)(z + 0.5f);
        float zf = z + 0.5f - zi;
        float zd0 = voxOutputArray[x][y][zi - 1];
        float zd1 = voxOutputArray[x][y][zi];
        float zd2 = voxOutputArray[x][y][zi + 1];
        res[2] = (zd1 - zd0) * (1.0f - zf) + (zd2 - zd1) * zf; // lerp
        return res;
    }

    // Overwrite floor voxels to a fixed pressure
    private void overwriteFloor(int floorPos){
        for(int x = 0; x < size; x++){
            for(int y = 0; y < size; y++){
                if(y > floorPos) continue;
                for(int z = 0; z < size; z++){
                    voxOutputArray[x][y][z] = groundPressure;
                }
            }
        }
    }

    // Overwrite air vocels to a fixed pressure
    private void overwriteAir(int airPos){
        for(int x = 0; x < size; x++){
            for(int y = 0; y < size; y++){
                for(int z = 0; z < size; z++){
                    /* Gladino dobimo tako, da Y dimenziji (višini)
                    * odštejemo prostor za valove (2 * airpos)
                    * ter prištejemo vrednost funkcije valovne površine,
                    *  ki se giblje od -airpos do airpos. */
                    int surface = this.size - 2 * airPos +
                            (int)Math.round(airPos/5 * (Math.sin(0.5 * x) + Math.cos(0.5 * z)));
                    if(y > surface)
                        voxOutputArray[x][y][z] = airPressure;
                }
            }
        }
    }

    // Append current vocOutbutArray to final output byte array
    private byte[] outputBytes;
    void appendToOutputFile(int imageNum, int currentIter){
        if(outputBytes == null)
            outputBytes = new byte[size*size*size*4*imageNum];

        int count = size * size * size * 4 * currentIter;

        for(int x = 0; x < size; x++){
            for(int y = 0; y < size; y++){
                for(int z = 0; z < size; z++){
                    outputBytes[count] = (byte)voxOutputArray[x][y][z];
                    count++;
                    // System.out.println(voxOutputArray[x][y][z]);
                    // System.out.println("Float val: " + voxOutputArray[x][y][z]);
                    if(true){
                        float[] grad = computeGradientsScript(x, y, z);
                        //float[] grad2 = computeGradients(x, y, z);
                        outputBytes[count] = (byte)grad[0];
                        outputBytes[count + 1] = (byte)grad[1];
                        outputBytes[count + 2] = (byte)grad[2];
                        count += 3;
                    }
                }
            }
        }
        //System.out.println(voxOutputArray[60][60][60]);
        System.out.println(".");
    }

    // Create RAW file from outputBytes array
    void createOutputFile(int id) throws IOException {
            FileOutputStream fout = new FileOutputStream(".\\generateVolume.out");
            fout.write(outputBytes);
            fout.close();
            System.out.println("Successfully generated output file.");
    }
}
