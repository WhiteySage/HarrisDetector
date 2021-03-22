import java.awt.*;
import java.awt.color.ColorSpace;
import java.awt.image.BufferedImage;
import java.awt.image.ColorConvertOp;
import java.io.IOException;
import java.util.ArrayList;
import java.awt.Color;
import java.util.Map;

public class Harris {


    /*
    Получаем чб изображение
     */
    public static BufferedImage convertToGrayScale(BufferedImage img) {
        ColorSpace cs = ColorSpace.getInstance(ColorSpace.CS_GRAY);
        ColorConvertOp op = new ColorConvertOp(cs, null);
        BufferedImage image = op.filter(img, null);
        return image;
    }

    /*
    Преобразовываем изображение в в массив из 0 и 1 для большей точности
     */
    public static double[][] getArrayFromBinaryImage(BufferedImage bufferedImage) {
        BufferedImage source;
        double[][] resultImage = new double[0][0];
        try {
            source = bufferedImage;
            int rowCount = source.getHeight();
            int colCount = source.getWidth();
            resultImage = new double[rowCount][colCount];
            for (int i = 0; i < rowCount; i++) {
                for (int j = 0; j < colCount; j++) {
                    Color pixelColor = new Color(source.getRGB(j, i));
                    int red = (int) (pixelColor.getRed());
                    int green = (int) (pixelColor.getGreen());
                    int blue = (int) (pixelColor.getBlue());

                    /*
                    240 и выше - это белый
                     */
                    if (red > 240 && green > 240 && blue > 240)
                        resultImage[i][j] = 1;
                    else
                        resultImage[i][j] = 0;
                }
            }

        } catch (Exception e) {
            e.printStackTrace();
        }


        return resultImage;
    }


    /*
    Выводим матрицу на экран
    */
    public static void printMatrix(double[][] matrix) {
        String result = "";
        for (int i = 0; i < getRowCount(matrix); i++) {
            for (int j = 0; j < getColumnCount(matrix); j++) {
                result += matrix[i][j] + " ";
            }
            result += "\n";
        }
        System.out.println(result);
    }

    /*
    Возвращает количество столбцов
     */
    private static int getColumnCount(double[][] matrix) {
        return matrix[0].length;
    }

    /*
    Возвращает количество строк
    */
    private static int getRowCount(double[][] matrix) {
        return matrix.length;
    }


    /*
    Создание ядра
     */
    public static double[][] createnKernel(double sigma) {
        int W = 3;
        double[][] kernel = new double[W][W];

        double mean = W / 2;
        double sum = 0;

        for (int x = 0; x < W; ++x) {
            for (int y = 0; y < W; ++y) {

                kernel[x][y] = (Math.exp(-0.5 * (Math.pow((x - mean) / sigma, 2.0) + Math.pow((y - mean) / sigma, 2.0)))
                        / (2 * Math.PI * sigma * sigma));

                sum += kernel[x][y];
            }
        }

        for (int i = 0; i < W; i++) {
            for (int j = 0; j < W; j++) {
                kernel[i][j] /= sum;
            }
        }

        return kernel;
    }

    /*
    Получаем DoG
     */
    public static double[][] GausX(double sigma) {
        double[] filter = {-1, 0, 1};
        double[][] kernel = createnKernel(sigma);
        for (int i = 0; i < kernel.length; i++) {
            for (int j = 0; j < kernel.length; j++) {
                kernel[i][j] *= filter[j];
            }
        }
        return kernel;
    }

    public static double[][] GausY(double sigma) {
        double[] filter = {-1, 0, 1};
        double[][] kernel = createnKernel(sigma);
        for (int i = 0; i < kernel.length; i++) {
            for (int j = 0; j < kernel.length; j++) {
                kernel[i][j] *= filter[i];
            }
        }
        return kernel;
    }

    /*
    Свертка, Ix, Iy
     */
    public static double[][] Convolution(double[][] original, double[][] kernel) {

        int additional = kernel.length - 1;
        int semiadittional = additional / 2;

        double[][] newMatrix = new double[original[0].length + additional][original.length + additional];


        //Big matrix
        for (int x = 0, i = semiadittional + x; x < original[0].length; x++, i++) {
            for (int y = 0, j = semiadittional + y; y < original.length; y++, j++) {
                newMatrix[i][j] =
                        original[y][x];
            }
        }

        double[][] result = new double[original.length][original[0].length];

        for (int x = 0, i = semiadittional + x; x < original[0].length; x++, i++) {
            for (int y = 0, j = semiadittional + y; y < original.length; y++, j++) {
                double[][] sub = SubMatrix(i - semiadittional, j - semiadittional, i + semiadittional, j + semiadittional, newMatrix);
                result[y][x] = ConvolutionElements(sub, kernel);
            }
        }


        return result;
    }

    public static double[][] SubMatrix(int x0, int y0, int xN, int yN, double[][] original) {
        double[][] result = new double[xN - x0 + 1][yN - y0 + 1];

        for (int x = 0, i = x0 + x; x < result[0].length; x++, i++) {
            for (int y = 0, j = y0 + y; y < result.length; y++, j++) {
                result[x][y] = original[i][j];
            }
        }
        return result;
    }

    public static int ConvolutionElements(double[][] matrix, double[][] kernel) {

        int result = 0;
        for (int i = 0; i < kernel[0].length; i++) {
            for (int j = 0; j < kernel.length; j++) {
                result += matrix[i][j] * kernel[i][j];
            }
        }
        return result;

    }


    /*
   Получение изображения из массива
    */
    public static BufferedImage drawImage(double[][] matrix) {

        BufferedImage img = new BufferedImage(matrix.length, matrix[0].length, BufferedImage.TYPE_INT_RGB);

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                int pixel = (int) matrix[i][j];
                img.setRGB(i, j, pixel);
            }
        }

        return img;
    }

    public static int getGrayScale(int rgb) {
        int r = (rgb >> 16) & 0xff;
        int g = (rgb >> 8) & 0xff;
        int b = (rgb) & 0xff;


        int gray = (int) (0.2126 * r + 0.7152 * g + 0.0722 * b);

        return gray;
    }


    /*/
    Преобразовение изображения в массив
     */
    static double[][] transformImageToArray(BufferedImage bufferedImage) {
        int width = bufferedImage.getWidth();
        int height = bufferedImage.getHeight();

        double[][] image = new double[height][width];

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                Color color = new Color(bufferedImage.getRGB(j, i));
                image[i][j] = color.getRed();
            }
        }
        return image;
    }


    /*
    Нормализируем матрицу
     */
    public static double[][] normalize(double[][] matrix) {
        double maxR = 0;
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (maxR < matrix[i][j]) {
                    maxR = matrix[i][j];
                }
            }
        }

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                matrix[i][j] /= maxR;
            }
        }


        return matrix;
    }





    /*
    Детектор углов Харриса
     */

    public static double[][] findR(BufferedImage image, double k) throws IOException {

        image = greyImage(image);

        //  Gx
        double[][] kernelX = Harris.GausX(0.5);
        //  Gy
        double[][] kernelY = Harris.GausY(0.5);

        //  Ix
        double[][] I_x = Harris.Convolution(Harris.transformImageToArray(image), kernelY);
        //  Iy
        double[][] I_y = Harris.Convolution(Harris.transformImageToArray(image), kernelX);

        // Ix^2
        double[][] I_x2 = new double[I_x.length][I_x[0].length];
        for (int i = 0; i < I_x.length; i++) {
            for (int j = 0; j < I_x[0].length; j++) {
                I_x2[i][j] = I_x[i][j] * I_x[i][j];
            }
        }
        // Iy^2
        double[][] I_y2 = new double[I_y.length][I_y[0].length];
        for (int i = 0; i < I_x.length; i++) {
            for (int j = 0; j < I_x[0].length; j++) {
                I_y2[i][j] = I_y[i][j] * I_y[i][j];
            }
        }
        // Ixy
        double[][] I_xy = new double[I_y.length][I_y[0].length];
        for (int i = 0; i < I_x.length; i++) {
            for (int j = 0; j < I_x[0].length; j++) {
                I_xy[i][j] = I_x[i][j] * I_y[i][j];
            }
        }
        // Фильтр с другой сигмой
        double[][] kernel = Harris.createnKernel(1);

        //Sx2
        double[][] S_x2 = Harris.Convolution(I_x2, kernel);


        //Sy2
        double[][] S_y2 = Harris.Convolution(I_y2, kernel);


        //Sxy
        double[][] S_xy = Harris.Convolution(I_xy, kernel);


        double low = 0.5;

        /*
        Матрицу R создается
         */
        double[][] R = new double[image.getHeight()][image.getWidth()];

        double maxR = 0;

        for (int i = 0; i < image.getHeight(); i++) {
            for (int j = 0; j < image.getWidth(); j++) {
                double trace = S_x2[i][j] + S_y2[i][j];
                double det = S_x2[i][j] * S_y2[i][j] - S_xy[i][j] * S_xy[i][j];
                double r = det - k * trace * trace;

                R[i][j] = r;
            }
        }

        return R;
    }


    //Зарисовка углов

//    public static Map<Integer, Integer[]> cornerDetector(double[][] R, BufferedImage image,int NonMaximaFactor) {
//
//        int nmf = NonMaximaFactor / 2;
//        int nmfm = nmf + 1;
//
//        HashSet<Integer[]> corners = new HashSet<>();
//        ArrayList<Integer> keys = new ArrayList<>();
//
//        for (int i = nmf; i < image.Rows - nmf; i++)
//        {
//            for (int j = nmf; j < GrayImage.Columns - nmf; j++)
//            {
//                float localMaximum = R.data[i, j];
//                if (R.data[i, j] > RFactor)
//                {
//                    int MaxC = j * 1, MaxR = i * 1;
//                    for (int x = -nmf; x < nmfm; x++)
//                    {
//                        for (int y = -nmf; y < nmfm; y++)
//                        {
//                            int currentR = i + x, currentC = j + y;
//                            float current = R.data[currentR, currentC];
//                            if (current > 0)
//                            {
//                                if (current > localMaximum)
//                                {
//                                    localMaximum = current;
//                                    MaxR = currentR;
//                                    MaxC = currentC;
//                                }
//                                else if (current < localMaximum)
//                                {
//                                    R.data[currentR, currentC] = 0;
//                                }
//                            }
//                        }
//                    }
//                    int key = (MaxR << 16) + MaxC;
//                    if (!keys.Contains(key))
//                    {
//                        corners.Add(new int[] { MaxR, MaxC });
//                        keys.Add(key);
//                    }
//                }
//            }
//        }
//        return (R, corners);
//
//    }


     /*
    Область вокруг интересных точек
     */
    /*

    public static ArrayList<Double[]> areaAroundInterestingPoints(double[][] R, Map<Integer, Integer[]> corners, BufferedImage image) {

        ArrayList<Double[]> arrayList = new ArrayList<>();

        for (Map.Entry<Integer, Integer[]> element : corners.entrySet()) {

            Integer[] coordinate = element.getValue();
            int i = coordinate[0];
            int j = coordinate[1];
            Double[] current = new Double[25];
            int index = 0;

            for (int y = i - 2; y < i + 3; y++) {
                for (int x = j - 2; x < j + 3; x++, index++) {

                    if (x > 0 && y > 0 && x < image.getHeight() && y < image.getWidth()) {
                        current[index] = R[x][y];

                    } else {
                        current[index] = 0.0;
                    }

                }
            }

            arrayList.add(current);

        }
        return arrayList;
    }


    */


    /*
    Преобразуем в диапазон 0-1
     */
    public static double[][] threshhold(double[][] R) {
        double thresh = 0.001;
        for (int i = 0; i < R.length; i++) {
            for (int j = 0; j < R[0].length; j++) {
                if (R[i][j] < thresh) {
                    R[i][j] = 0;
                }
            }
        }
        return R;
    }

    /*
    Область вокруг интересных точек
    Находим угол и в области угла 5х5 находим  максимальное делаем и присваеваем остальным 0
     */
    public static double[][] cornerDetector2(double[][] R, BufferedImage image) {
        double max;
        int Y;
        int X;
        for (int i = 2; i < R.length - 2; i++) {
            for (int j = 2; j < R[0].length - 2; j++) {
                max = 0;
                Y = 0;
                X = 0;
                if (R[i][j] > 0) {
                    for (int y = i - 2; y < i + 3; y++) {
                        for (int x = j - 2; x < j + 3; x++) {
                            if (max < R[y][x]) {
                                max = R[y][x];
                                Y = y;
                                X = x;
                                R[y][x] = 0;
                            } else {
                                R[y][x] = 0;
                            }
                        }
                    }
                    R[Y][X] = 1;
                    Color c = new Color(255, 0, 0);
                    image.setRGB(X, Y, c.getRGB());
                }
            }
        }

        return R;
    }


//    public static ArrayList<Integer[]> areaAroundInterstingPoints2(double[][] R, BufferedImage image) {
//
//        ArrayList<Integer[]> list = new ArrayList<>();
//        Integer[] array = new Integer[25];
//        int index = 0;
//        for (int i = 2; i < R.length - 2; i++) {
//            for (int j = 2; j < R[0].length - 2; j++) {
//                if (R[i][j] > 0) {
//                    for (int y = i - 2; y < i + 3; y++) {
//                        for (int x = j - 2; x < j + 3; x++) {
//
//                            Color c = new Color(image.getRGB(x, y));
//                            Integer rgb = c.getRed();
//                            array[index] = rgb;
//                            index++;
//
//
//                        }
//                    }
//                    list.add(array);
//                    index = 0;
//                }
//            }
//        }
//        return list;
//    }

    /*
    Преобразуем изображение в серый
     */
    public static BufferedImage greyImage(BufferedImage image) throws IOException {
        try {

            int width = image.getWidth();
            int height = image.getHeight();
            for (int i = 0; i < width; i++) {
                for (int j = 0; j < height; j++) {
                    Color c = new Color(image.getRGB(i, j));

                    int red = (int) (c.getRed() * 0.2126);
                    int green = (int) (c.getGreen() * 0.7152);
                    int blue = (int) (c.getBlue() * 0.0722);

                    Color newColor = new Color(red + green + blue,
                            red + green + blue, red + green + blue);
                    image.setRGB(i, j, newColor.getRGB());
                }
            }
        } catch (Exception e) {
            e.getMessage();
        }
        return image;
    }


    /*
    Вывод матрицы дискриптора
     */
    public static void printPatch(ArrayList<Double[]> descriptor) {
        for (int i = 0; i < descriptor.size(); i++) {
            for (int z = 0; z < descriptor.get(i).length; z++) {
                System.out.print(descriptor.get(i)[z] + " ");
            }
            System.out.println();
        }
    }

    /*
    Вывод матрицы углов
    */
    public static void printCorners(ArrayList<Integer[]> descriptor) {
        for (int i = 0; i < descriptor.size(); i++) {
            for (int z = 0; z < descriptor.get(i).length; z++) {
                System.out.print(descriptor.get(i)[z] + " ");
            }
            System.out.println();
        }
    }


    /*
    Заносим координаты углов в массив
     */
    public static ArrayList<Integer[]> findCorners(double[][] R) {
        ArrayList<Integer[]> corners = new ArrayList<>();
        for (int i = 0; i < R.length; i++) {
            for (int j = 0; j < R[0].length; j++) {
                if (R[i][j] > 0) {
                    corners.add(new Integer[]{i, j});
                }
            }
        }
        return corners;
    }

    /*
    Область яркости пикселей вокруг интересных точек для дескриптора
     */
    public static ArrayList<Double[]> areaAroundInterstingPoints2(ArrayList<Integer[]> corners, BufferedImage image, int descriptorFactor) {

        ArrayList<Double[]> descriptors = new ArrayList<Double[]>();
        int descriptorFactor2 = descriptorFactor * descriptorFactor;
        int semiFactor = descriptorFactor / 2;
        int semiFactorIncr = semiFactor + 1;
        semiFactor = 0 - semiFactor;

        for (Integer[] corner : corners) {
            int r = corner[0], c = corner[1];
            int index = 0;
            Double[] descriptor = new Double[descriptorFactor2];
            for (int x = semiFactor; x < semiFactorIncr; x++) {
                for (int y = semiFactor; y < semiFactorIncr; y++, index++) {
                    int currentR = r + x;
                    int currentC = c + y;
                    if (currentR > 0 && currentR < image.getWidth() && currentC > 0 && currentC < image.getHeight()) {
                        Color color = new Color(image.getRGB(currentR, currentC));
                        descriptor[index] = Double.valueOf(color.getRed());
                    } else {
                        descriptor[index] = 0.0;
                    }
                }
            }
            descriptors.add(descriptor);
        }
        return descriptors;

    }

    /*
    Сравниваем дескрипторы
     */
    public static ArrayList<Integer[]> matchDescriptor(ArrayList<Double[]> descriptors1, ArrayList<Double[]> descriptors2, int descriptorFactor) {

        Double[][] comparison = new Double[descriptors1.size()][descriptors2.size()];
        ArrayList<Integer[]> pairs = new ArrayList<>();
        int descriptorFactor2 = descriptorFactor * descriptorFactor;

        for (int r = 0; r < comparison.length; r++) {
            for (int c = 0; c < comparison[0].length; c++) {
                float sqrt = 0;
                for (int index = 0; index < descriptorFactor2; index++) {
                    Double delta = descriptors1.get(r)[index] - descriptors2.get(c)[index];
                    sqrt += delta * delta;
                }
                comparison[r][c] = Math.sqrt(sqrt);
            }
        }

        for (int r = 0; r < comparison.length; r++) {
            Double[] currentMin = new Double[]{Double.MAX_VALUE, Double.MAX_VALUE};
            int minId = 0;
            for (int c = 0; c < comparison[0].length; c++) {
                if (comparison[r][c] < currentMin[0]) {
                    currentMin[1] = currentMin[0];
                    currentMin[0] = comparison[r][c];
                    minId = c;
                }
            }
            if (currentMin[0] / currentMin[1] < 0.8f) {
                for (int i = 0; i < comparison[0].length; i++) {
                    comparison[i][minId] = Double.POSITIVE_INFINITY;
                }
                pairs.add(new Integer[]{r, minId});
            }
        }


        return pairs;
    }

    /*
    Рисуем линию между углами, которые совпали
     */
    public static BufferedImage drawMatches(BufferedImage image1, BufferedImage image2, ArrayList<Integer[]> pairs, ArrayList<Integer[]> coordinate1, ArrayList<Integer[]> coordinate2) {

        BufferedImage image = new BufferedImage(image1.getWidth() + image2.getWidth(), Math.max(image1.getHeight(), image2.getHeight()), BufferedImage.TYPE_INT_RGB);


        for (int x = 0; x < image1.getWidth(); x++) {
            for (int y = 0; y < image1.getHeight(); y++) {
                image.setRGB(x, y, image1.getRGB(x, y));
            }
        }

        for (int x = image2.getWidth(), x0 = 0; x0 < image2.getWidth(); x++, x0++) {
            for (int y = 0; y < image2.getHeight(); y++) {
                image.setRGB(x, y, image2.getRGB(x0, y));
            }
        }

        Graphics2D g2d = image.createGraphics();
        Color color = new Color(0, 255, 0);


        for (Integer[] pair : pairs) {
            Integer[] c1 = coordinate1.get(pair[0]);
            Integer[] c2 = coordinate2.get(pair[1]);

            Point p1 = new Point(c1[1], c1[0]);
            Point p2 = new Point(c2[1] + image1.getWidth(), c2[0]);
            g2d.setColor(color);
            g2d.drawLine(p1.x, p1.y, p2.x, p2.y);

            image.setRGB(p1.x, p1.y, color.getRGB());
            image.setRGB(p2.x, p2.y, color.getRGB());
//            g2d.fillOval(p1.x-2,p1.y-2,5,5);
//            g2d.fillOval(p2.x-2,p2.y-2,5,5);
        }

        return image;
    }


}















