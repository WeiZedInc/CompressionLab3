namespace CompressionLab3
{
    public partial class MainFormCompressor : Form
    {
        #region Windows Form Designer generated code

        private System.ComponentModel.IContainer components = null;
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }
        private void InitializeComponent()
        {
            this.components = new System.ComponentModel.Container();
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(800, 450);
            this.Text = "Form1";
        }

        #endregion

        private Bitmap originalImage;
        private Bitmap dctImage;
        private Bitmap dwtImage;

        public MainFormCompressor(EnvironmentVariableTarget a)
        {
            InitializeComponent();

            //LoadImage();
            CreateTestImage();

            ApplyDCT();
            ApplyDWT();
            ShowImages();
            DisplayImageSizes();
        }

        private void LoadImage() => originalImage = new Bitmap("path_to_your_image.jpg");
        private void CreateTestImage()
        {
            int size = 512;
            originalImage = new Bitmap(size, size);
            using (Graphics g = Graphics.FromImage(originalImage))
            {
                g.Clear(Color.White); // Білий фон

                // Чорний квадрат у верхньому лівому куті
                g.FillRectangle(Brushes.Black, 50, 50, 150, 150);

                // Сірий квадрат у правому нижньому куті
                g.FillRectangle(Brushes.Gray, 300, 300, 150, 150);

                // Світло-сіре коло в центрі
                g.FillEllipse(Brushes.LightGray, 200, 200, 100, 100);
            }
        }

        private void ApplyDCT()
        {
            int blockSize = 8;
            int width = originalImage.Width;
            int height = originalImage.Height;
            dctImage = new Bitmap(width, height);

            for (int x = 0; x < width; x += blockSize)
            {
                for (int y = 0; y < height; y += blockSize)
                {
                    double[,] block = GetBlock(originalImage, x, y, blockSize);
                    double[,] dctBlock = DCTTransform(block);
                    dctBlock = Quantize(dctBlock);
                    double[,] inverseDCTBlock = InverseDCTTransform(dctBlock);
                    SetBlock(dctImage, inverseDCTBlock, x, y);
                }
            }
        }
        private void ApplyDWT()
        {
            int width = originalImage.Width;
            int height = originalImage.Height;
            dwtImage = new Bitmap(width, height);
            double[,] imageMatrix = ImageToMatrix(originalImage);
            double[,] dwtMatrix = HaarWaveletTransform(imageMatrix);
            dwtMatrix = Quantize(dwtMatrix);
            dwtImage = MatrixToImage(dwtMatrix);
        }

        private double[,] DCTTransform(double[,] block)
        {
            int N = block.GetLength(0);
            double[,] dct = new double[N, N];
            for (int u = 0; u < N; u++)
            {
                for (int v = 0; v < N; v++)
                {
                    double sum = 0.0;
                    for (int x = 0; x < N; x++)
                    {
                        for (int y = 0; y < N; y++)
                            sum += block[x, y] * Math.Cos((2 * x + 1) * u * Math.PI / (2 * N)) * Math.Cos((2 * y + 1) * v * Math.PI / (2 * N));
                    }
                    dct[u, v] = sum * (u == 0 ? 1 / Math.Sqrt(2) : 1) * (v == 0 ? 1 / Math.Sqrt(2) : 1);
                }
            }
            return dct;
        }
        private double[,] InverseDCTTransform(double[,] dctBlock)
        {
            int N = dctBlock.GetLength(0);
            double[,] block = new double[N, N];
            for (int x = 0; x < N; x++)
            {
                for (int y = 0; y < N; y++)
                {
                    double sum = 0.0;
                    for (int u = 0; u < N; u++)
                    {
                        for (int v = 0; v < N; v++)
                        {
                            sum += (u == 0 ? 1 / Math.Sqrt(2) : 1) * (v == 0 ? 1 / Math.Sqrt(2) : 1) *
                                   dctBlock[u, v] * Math.Cos((2 * x + 1) * u * Math.PI / (2 * N)) * Math.Cos((2 * y + 1) * v * Math.PI / (2 * N));
                        }
                    }
                    block[x, y] = Math.Round(sum / 4); // нормалізація
                }
            }
            return block;
        }

        private double[,] HaarWaveletTransform(double[,] image)
        {
            int N = image.GetLength(0);
            double[,] transformed = new double[N, N];

            for (int i = 0; i < N; i += 2)
            {
                for (int j = 0; j < N; j += 2)
                {
                    double avg = (image[i, j] + image[i + 1, j] + image[i, j + 1] + image[i + 1, j + 1]) / 4;
                    transformed[i / 2, j / 2] = avg;
                    transformed[i / 2, j / 2 + N / 2] = image[i, j] - avg;
                    transformed[i / 2 + N / 2, j / 2] = image[i + 1, j] - avg;
                    transformed[i / 2 + N / 2, j / 2 + N / 2] = image[i, j + 1] - avg;
                }
            }
            return transformed;
        }
        private double[,] Quantize(double[,] matrix)
        {
            int N = matrix.GetLength(0);
            double quantizationFactor = 10.0; // агресивність
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (i > N / 2 || j > N / 2)
                        matrix[i, j] = 0; // Обнулення високочастотних компонентів
                    else
                        matrix[i, j] = Math.Round(matrix[i, j] / quantizationFactor) * quantizationFactor;
                }
            }
            return matrix;
        }

        private double CalculatePSNR(Bitmap original, Bitmap compressed)
        {
            double mse = 0;
            int width = original.Width;
            int height = original.Height;
            for (int x = 0; x < width; x++)
            {
                for (int y = 0; y < height; y++)
                {
                    double diff = original.GetPixel(x, y).R - compressed.GetPixel(x, y).R;
                    mse += diff * diff;
                }
            }
            mse /= (width * height);
            return 10 * Math.Log10(255 * 255 / mse);
        }

        private void ShowImages()
        {
            int imageSize = 256;
            int gap = 20;

            Label originalLabel = new Label
            {
                Text = "Оригінал",
                AutoSize = true,
                Location = new Point(10 + (imageSize - TextRenderer.MeasureText("Оригінал", SystemFonts.DefaultFont).Width) / 2, 10)
            };
            PictureBox originalBox = new PictureBox
            {
                Image = originalImage,
                Size = new Size(imageSize, imageSize),
                Location = new Point(10, originalLabel.Bottom + 5),
                BorderStyle = BorderStyle.FixedSingle,
                SizeMode = PictureBoxSizeMode.Zoom
            };

            Label dctLabel = new Label
            {
                Text = "DCT",
                AutoSize = true,
                Location = new Point(originalBox.Right + gap + (imageSize - TextRenderer.MeasureText("DCT", SystemFonts.DefaultFont).Width) / 2, 10)
            };
            PictureBox dctBox = new PictureBox
            {
                Image = dctImage,
                Size = new Size(imageSize, imageSize),
                Location = new Point(originalBox.Right + gap, dctLabel.Bottom + 5),
                BorderStyle = BorderStyle.FixedSingle,
                SizeMode = PictureBoxSizeMode.Zoom
            };

            Label dwtLabel = new Label
            {
                Text = "DWT",
                AutoSize = true,
                Location = new Point(dctBox.Right + gap + (imageSize - TextRenderer.MeasureText("DWT", SystemFonts.DefaultFont).Width) / 2, 10)
            };
            PictureBox dwtBox = new PictureBox
            {
                Image = dwtImage,
                Size = new Size(imageSize, imageSize),
                Location = new Point(dctBox.Right + gap, dwtLabel.Bottom + 5),
                BorderStyle = BorderStyle.FixedSingle
            };

            this.Controls.Add(originalLabel);
            this.Controls.Add(originalBox);
            this.Controls.Add(dctLabel);
            this.Controls.Add(dctBox);
            this.Controls.Add(dwtLabel);
            this.Controls.Add(dwtBox);

            // Обчислення PSNR
            double psnrDCT = CalculatePSNR(originalImage, dctImage);
            double psnrDWT = CalculatePSNR(originalImage, dwtImage);

            Label psnrDCTLabel = new Label
            {
                Text = $"PSNR для DCT: {psnrDCT:F2} дБ",
                AutoSize = true,
                Location = new Point(dctBox.Left + (imageSize - TextRenderer.MeasureText($"PSNR для DCT: {psnrDCT:F2} дБ", SystemFonts.DefaultFont).Width) / 2, dctBox.Bottom + 5)
            };

            Label psnrDWTLabel = new Label
            {
                Text = $"PSNR для DWT: {psnrDWT:F2} дБ",
                AutoSize = true,
                Location = new Point(dwtBox.Left + (imageSize - TextRenderer.MeasureText($"PSNR для DWT: {psnrDWT:F2} дБ", SystemFonts.DefaultFont).Width) / 2, dwtBox.Bottom + 5)
            };

            this.Controls.Add(psnrDWTLabel);
            this.Controls.Add(psnrDCTLabel);
        }
        private void DisplayImageSizes()
        {
            long originalSize = GetImageSizeInBytes(originalImage);
            long dctSize = GetImageSizeInBytes(dctImage);
            long dwtSize = GetImageSizeInBytes(dwtImage);

            double dctCompressionPercentage = 100.0 * (originalSize - dctSize) / originalSize;
            double dwtCompressionPercentage = 100.0 * (originalSize - dwtSize) / originalSize;

            int originalWidth = originalImage.Width;
            int originalHeight = originalImage.Height;

            int dctWidth = dctImage.Width;
            int dctHeight = dctImage.Height;

            int dwtWidth = dwtImage.Width;
            int dwtHeight = dwtImage.Height;

            Label originalSizeLabel = new Label
            {
                Text = $"Оригінал: {originalWidth}x{originalHeight}, {originalSize / 1024.0:F2} KB",
                AutoSize = true,
                Location = new Point(10, 400)
            };
            this.Controls.Add(originalSizeLabel);

            Label dctSizeLabel = new Label
            {
                Text = $"Після DCT: {dctWidth}x{dctHeight}, {dctSize / 1024.0:F2} KB ({dctCompressionPercentage:F2}%)",
                AutoSize = true,
                Location = new Point(270, 400)
            };
            this.Controls.Add(dctSizeLabel);

            Label dwtSizeLabel = new Label
            {
                Text = $"Після DWT: {dwtWidth}x{dwtHeight}, {dwtSize / 1024.0:F2} KB ({dwtCompressionPercentage:F2}%)",
                AutoSize = true,
                Location = new Point(530, 400)
            };
            this.Controls.Add(dwtSizeLabel);
        }

        private long GetImageSizeInBytes(Bitmap image, long quality = 50L)
        {
            using (MemoryStream ms = new MemoryStream())
            {
                var encoderParameters = new System.Drawing.Imaging.EncoderParameters(1);
                encoderParameters.Param[0] = new System.Drawing.Imaging.EncoderParameter(
                    System.Drawing.Imaging.Encoder.Quality, quality);

                var jpegEncoder = GetEncoder(System.Drawing.Imaging.ImageFormat.Jpeg);
                image.Save(ms, jpegEncoder, encoderParameters);
                return ms.Length;
            }
        }
        private System.Drawing.Imaging.ImageCodecInfo GetEncoder(System.Drawing.Imaging.ImageFormat format)
        {
            var codecs = System.Drawing.Imaging.ImageCodecInfo.GetImageDecoders();
            foreach (var codec in codecs)
            {
                if (codec.FormatID == format.Guid)
                    return codec;
            }
            return null;
        }

        private double[,] GetBlock(Bitmap image, int startX, int startY, int blockSize)
        {
            double[,] block = new double[blockSize, blockSize];
            for (int x = 0; x < blockSize; x++)
            {
                for (int y = 0; y < blockSize; y++)
                    block[x, y] = image.GetPixel(startX + x, startY + y).R;
            }
            return block;
        }
        private void SetBlock(Bitmap image, double[,] block, int startX, int startY)
        {
            int blockSize = block.GetLength(0);
            for (int x = 0; x < blockSize; x++)
            {
                for (int y = 0; y < blockSize; y++)
                {
                    int value = (int)Math.Clamp(block[x, y], 0, 255);
                    image.SetPixel(startX + x, startY + y, Color.FromArgb(value, value, value));
                }
            }
        }

        private double[,] ImageToMatrix(Bitmap image)
        {
            int width = image.Width;
            int height = image.Height;
            double[,] matrix = new double[width, height];
            for (int x = 0; x < width; x++)
            {
                for (int y = 0; y < height; y++)
                    matrix[x, y] = image.GetPixel(x, y).R;
            }
            return matrix;
        }
        private Bitmap MatrixToImage(double[,] matrix)
        {
            int width = matrix.GetLength(0);
            int height = matrix.GetLength(1);
            Bitmap image = new Bitmap(width, height);
            for (int x = 0; x < width; x++)
            {
                for (int y = 0; y < height; y++)
                {
                    int value = (int)Math.Clamp(matrix[x, y], 0, 255);
                    image.SetPixel(x, y, Color.FromArgb(value, value, value));
                }
            }
            return image;
        }
    }
}
