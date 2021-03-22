import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;

public class Main {
    public static void main(String[] args) throws IOException {

        // Коэфициент Харриса
        double k = 0.1;
        // Вот тут это область для дискриптора вокруг угла
        int factor = 55;

        // Коллекция для дескрептира
        ArrayList<Double[]> descriptorA = new ArrayList<>();
        ArrayList<Double[]> descriptorB = new ArrayList<>();

        // Чтение файла
        File file = new File("test1.jpg");
        File file1 = new File("test2.jpg");
        BufferedImage image = ImageIO.read(file);
        BufferedImage image1 = ImageIO.read(file1);


        // Найдем значение R для изображения 1
        double[][] R = Harris.findR(image, k);
        R = Harris.normalize(R);
        R = Harris.threshhold(R);

        // Найдем значение R для изображения 2
        double[][] R1 = Harris.findR(image1, k);
        R1 = Harris.normalize(R1);
        R1 = Harris.threshhold(R1);


        // Найдем углы на изображении
        R = Harris.cornerDetector2(R, image);
        R1 = Harris.cornerDetector2(R1, image1);


        // Найдем координаты углов
        var cornerA = Harris.findCorners(R);
        var cornerB = Harris.findCorners(R1);


        // Патч-дескриптор
        descriptorA = Harris.areaAroundInterstingPoints2(cornerA, image, factor);
        descriptorB = Harris.areaAroundInterstingPoints2(cornerB, image1, factor);

        // Найдем совпадения на изображениях
        var matches = Harris.matchDescriptor(descriptorA, descriptorB, factor);

        // Сопоставление точек интересов 2-х изображений на 1
        BufferedImage finalImage = Harris.drawMatches(image, image1, matches, cornerA, cornerB);
        File finalImg = new File("finalImg.png");
        ImageIO.write(finalImage, "png", finalImg);


        System.out.println("DONE");
    }


}
