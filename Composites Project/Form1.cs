using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Composites_Project
{
    public partial class Form1 : Form
    {
        List<compositeLayer> layers = new List<compositeLayer>();
        List<TextBox> angles = new List<TextBox>();
        int layerToDisplay = 0;
        Matrix<double> laminate_strains;
        Matrix<double> laminate_curvatures;


        public Form1()
        {
            InitializeComponent();
            angles.Add(layer1Angle_textBox);
            angles.Add(layer2Angle_textBox);
            angles.Add(layer3Angle_textBox);
            angles.Add(layer4Angle_textBox);
            angles.Add(layer5Angle_textBox);
        }

        private void textBoxTextChanged(object sender, EventArgs e)
        {
            layers.Clear();
            try
            {
                /*Assign composite material properties as static members*/
                compositeLayer.E1 = Convert.ToDouble(E1_textBox.Text);
                compositeLayer.E2 = Convert.ToDouble(E2_textBox.Text);
                compositeLayer.G12 = Convert.ToDouble(G12_textBox.Text);
                compositeLayer.v12 = Convert.ToDouble(v12_textBox.Text);
                compositeLayer.numberOfLayers = Convert.ToInt16(numLayers_textBox.Text);
                compositeLayer.layerThickness = Convert.ToDouble(plyThickness_textBox.Text);
                
                /*Assign the fores and moments as static members to compositeLayer.*/
                compositeLayer.Nx = Convert.ToDouble(Nx_textBox.Text);
                compositeLayer.Ny = Convert.ToDouble(Ny_textBox.Text);
                compositeLayer.Nxy = Convert.ToDouble(Nxy_textBox.Text);
                compositeLayer.Mx = Convert.ToDouble(Mx_textBox.Text);
                compositeLayer.My = Convert.ToDouble(My_textBox.Text);
                compositeLayer.Mxy = Convert.ToDouble(Mxy_textBox.Text);

                /* Allocate new compositeLayer classes for each layer in the part. */
                for (int i = 0; i < compositeLayer.numberOfLayers; i++)
                {
                    compositeLayer newLayer = new compositeLayer();
                    layers.Add(newLayer);
                    compositeLayer.midplaneLayer = i;
                }

                /* Assign an angle to each layer in the layers list. */
                for (int i = 0; i < layers.Count; i++)
                {
                    try
                    {
                        layers[i].t = Convert.ToDouble(angles[i].Text) * (Math.PI / 180);
                    }
                    catch { }

                    /* Determine the top and bottom of each layer. */
                    layers[i].top = (layers.Count - i) * compositeLayer.layerThickness;
                    layers[i].bottom = (layers.Count - (i + 1)) * compositeLayer.layerThickness;
                }
            }
            catch
            {
                return;
            }

            /*A,B,D matrix variables.*/
            double q11a = 0;
            double q11b = 0;
            double q11d = 0;
            double q12a = 0;
            double q12b = 0;
            double q12d = 0;
            double q16a = 0;
            double q16b = 0;
            double q16d = 0;
            double q22a = 0;
            double q22b = 0;
            double q22d = 0;
            double q26a = 0;
            double q26b = 0;
            double q26d = 0;
            double q66a = 0;
            double q66b = 0;
            double q66d = 0;
            
            /*Calculate the A,B,D matrix values.*/
            for (int i = 0; i < layers.Count; i++)
            {
                q11a += layers[i].q11t() * layers[i].hk1();
                q11b += layers[i].q11t() * (layers[i].hk2() + -layers[i].hk2());
                q11d += layers[i].q11t() * (layers[i].hk3() * 2);

                q12a += layers[i].q12t() * layers[i].hk1();
                q12b += layers[i].q12t() * (layers[i].hk2() + -layers[i].hk2());
                q12d += layers[i].q12t() * (layers[i].hk3() * 2);

                q16a += layers[i].q16t() * layers[i].hk1();
                q16b += layers[i].q16t() * (layers[i].hk2() + -layers[i].hk2());
                q16d += layers[i].q16t() * (layers[i].hk3() * 2);

                q26a += layers[i].q26t() * layers[i].hk1();
                q26b += layers[i].q26t() * (layers[i].hk2() + -layers[i].hk2());
                q26d += layers[i].q26t() * (layers[i].hk3() * 2);

                q22a += layers[i].q22t() * layers[i].hk1();
                q22b += layers[i].q22t() * (layers[i].hk2() + -layers[i].hk2());
                q22d += layers[i].q22t() * (layers[i].hk3() * 2);

                q66a += layers[i].q66t() * layers[i].hk1();
                q66b += layers[i].q66t() * (layers[i].hk2() + -layers[i].hk2());
                q66d += layers[i].q66t() * (layers[i].hk3() * 2);
            }

            /*Assemble the ABD matrix(cies).*/
            Matrix<double> a_matrix = DenseMatrix.OfArray(new double[,]{ { 2 * q11a, 2 * q12a, 2 * q16a, },
                                                                        { 2 * q12a, 2 * q22a, 2 * q26a },
                                                                        { 2 * q16a, 2 * q26a, 2 * q66a } });
            Matrix<double> b_matrix = DenseMatrix.OfArray(new double[,]{ { q11b / 2, q12b / 2, q16b / 2, },
                                                                        { q12b / 2, q22b / 2, q26b / 2 },
                                                                        { q16b / 2, q26b / 2, q66b / 2 } });
            Matrix<double> d_matrix = DenseMatrix.OfArray(new double[,]{ { q11d / 3, q12d / 3, q16d / 3, },
                                                                        { q12d / 3, q22d / 3, q26d / 3 },
                                                                        { q16d / 3, q26d / 3, q66d / 3 } });
            Matrix<double> abd_matrix = DenseMatrix.OfArray(new double[,]{ { 2 * q11a, 2 * q12a, 2 * q16a, q11b / 2, q12b / 2, q16b / 2 },
                                                                            { 2 * q12a, 2 * q22a, 2 * q26a, q12b / 2, q22b / 2, q26b / 2 },
                                                                            { 2 * q16a, 2 * q26a, 2 * q66a, q16b / 2, q26b / 2, q66b / 2 },
                                                                            { q11b / 2, q12b / 2, q16b / 2, q11d / 3, q12d / 3, q16d / 3 },
                                                                            { q12b / 2, q22b / 2, q26b / 2, q12d / 3, q22d / 3, q26d / 3 },
                                                                            { q16b / 2, q26b / 2, q66b / 2, q16d / 3, q26d / 3, q66d / 3 } });
            Matrix<double> abd_inv_matrix = abd_matrix.Inverse();

            /*Assemble force and moment matrix.*/
            Matrix<double> forces_matrix = DenseMatrix.OfArray(new double[,] { { compositeLayer.Nx },
                                                                                { compositeLayer.Ny },
                                                                                { compositeLayer.Nxy },
                                                                                { compositeLayer.Mx },
                                                                                { compositeLayer.My },
                                                                                { compositeLayer.Mxy } });

            /*Solve for the strain and curvature matrix of the whole composite part.*/
            Matrix<double> stress_strain_matrix = abd_inv_matrix.Multiply(forces_matrix);
            
            /*Break the strain and curvature matrix apart into two seperate matricies, one for strain, one for curvatures.*/
            laminate_strains = DenseMatrix.OfArray(new double[,] { { stress_strain_matrix[0,0]},
                                                                    { stress_strain_matrix[1,0] },
                                                                    { stress_strain_matrix[2,0] }});
            laminate_curvatures = DenseMatrix.OfArray(new double[,] { { stress_strain_matrix[3, 0] },
                                                                    { stress_strain_matrix[4, 0] },
                                                                    { stress_strain_matrix[5, 0] }});

            for (int i = 0; i < compositeLayer.numberOfLayers; i++)
            {
                Matrix<double> qBarMatrix = DenseMatrix.OfArray(new double[,] { { layers[i].q11t(), layers[i].q12t(), layers[i].q16t()},
                                                                                { layers[i].q12t(), layers[i].q22t(), layers[i].q26t()},
                                                                                { layers[i].q16t(), layers[i].q26t(), layers[i].q66t()},});
                layers[i].top_stress_strain.strains = laminate_strains.Add(laminate_curvatures.Multiply(-layers[i].top));
                layers[i].bottom_stress_strain.strains = laminate_strains.Add(laminate_curvatures.Multiply(-layers[i].bottom));
                layers[i].top_stress_strain.stresses = qBarMatrix.Multiply(layers[i].top_stress_strain.strains);
                layers[i].bottom_stress_strain.stresses = qBarMatrix.Multiply(layers[i].bottom_stress_strain.strains);
                displayLayer(0);
            }
        }

        void displayLayer(int layerNumber)
        {
            try
            {
                compositeLayer thisLayer = layers[layerNumber];
                q11t_textBox.Text = Convert.ToString(thisLayer.q11t());
                q12t_textBox.Text = Convert.ToString(thisLayer.q12t());
                q12t_textBox2.Text = Convert.ToString(thisLayer.q12t());
                q16t_textBox.Text = Convert.ToString(thisLayer.q16t());
                q16t_textBox2.Text = Convert.ToString(thisLayer.q16t());
                q22t_textBox.Text = Convert.ToString(thisLayer.q22t());
                q26t_textBox.Text = Convert.ToString(thisLayer.q26t());
                q26t_textBox2.Text = Convert.ToString(thisLayer.q26t());
                q66t_textBox.Text = Convert.ToString(thisLayer.q66t());

                int roundingDigits = 6;
                exTop_textBox.Text = Convert.ToString(Math.Round(thisLayer.top_stress_strain.strains[0, 0], roundingDigits));
                eyTop_textBox.Text = Convert.ToString(Math.Round(thisLayer.top_stress_strain.strains[1, 0], roundingDigits));
                YxyTop_textBox.Text = Convert.ToString(Math.Round(thisLayer.top_stress_strain.strains[2, 0], roundingDigits));

                exBottom_textBox.Text = Convert.ToString(Math.Round(thisLayer.bottom_stress_strain.strains[0, 0], roundingDigits));
                eyBottom_textBox.Text = Convert.ToString(Math.Round(thisLayer.bottom_stress_strain.strains[1, 0], roundingDigits));
                YxyBottom_textBox.Text = Convert.ToString(Math.Round(thisLayer.bottom_stress_strain.strains[2, 0], roundingDigits));

                sigmaXtop_textBox.Text = Convert.ToString(Math.Round(thisLayer.top_stress_strain.stresses[0, 0], roundingDigits));
                sigmaYtop_textBox.Text = Convert.ToString(Math.Round(thisLayer.top_stress_strain.stresses[1, 0], roundingDigits));
                tauXYTop_textBox.Text = Convert.ToString(Math.Round(thisLayer.top_stress_strain.stresses[2, 0], roundingDigits));

                sigmaXBottom_textBox.Text = Convert.ToString(Math.Round(thisLayer.bottom_stress_strain.stresses[0, 0], roundingDigits));
                sigmaYBottom_textBox.Text = Convert.ToString(Math.Round(thisLayer.bottom_stress_strain.stresses[1, 0], roundingDigits));
                tauXYBottom_textBox.Text = Convert.ToString(Math.Round(thisLayer.bottom_stress_strain.stresses[2, 0], roundingDigits));
            }
            catch
            {
                MessageBox.Show("Please make sure all values are entered.");
            }

        }

        private void nextLayerButton_Click(object sender, EventArgs e)
        {
            if (layerToDisplay >= (compositeLayer.numberOfLayers - 1))
            {
                layerToDisplay = compositeLayer.numberOfLayers - 1;
            }
            else
            {
                layerToDisplay++;
            }
            displayedLayer_label.Text = "Layer " + Convert.ToString(layerToDisplay + 1);
            displayLayer(layerToDisplay);
        }

        private void previousLayerButton_Click(object sender, EventArgs e)
        {
            if (layerToDisplay == 0)
            {
                layerToDisplay = 0;
            }
            else
            {
                layerToDisplay--;
            }
            displayedLayer_label.Text = "Layer " + Convert.ToString(layerToDisplay + 1);
            displayLayer(layerToDisplay);
        }

        public static int[,] matrixMult(int[,] Matrix1, int[,] Matrix2)
        {
            if (Object.ReferenceEquals(null, Matrix1))
                throw new ArgumentNullException("Matrix1");
            else if (Object.ReferenceEquals(null, Matrix2))
                throw new ArgumentNullException("Matrix2");

            int r1 = Matrix1.GetLength(0);
            int c1 = Matrix1.GetLength(1);

            int r2 = Matrix2.GetLength(0);
            int c2 = Matrix2.GetLength(1);

            if (c1 != r2)
                throw new ArgumentOutOfRangeException("Matrix2", "Matrixes dimensions don't fit.");

            int[,] result = new int[r1, c2];

            // Naive matrix multiplication: O(n**3) 
            // Use Strassen algorithm O(n**2.81) in case of big matrices
            for (int r = 0; r < r1; ++r)
                for (int c = 0; c < c2; ++c)
                {
                    int s = 0;

                    for (int z = 0; z < c1; ++z)
                        s += Matrix1[r, z] * Matrix2[z, c];

                    result[r, c] = s;
                }

            return result;
        }

        private void Form1_Load(object sender, EventArgs e)
        {

        }
    }



    class compositeLayer
    {
        public static double E1;
        public static double E2;
        public static double G12;
        public static double v12;
        public static int numberOfLayers;
        public static double layerThickness;
        public static int midplaneLayer;
        public static double Nx;
        public static double Ny;
        public static double Nxy;
        public static double Mx;
        public static double My;
        public static double Mxy;
        static double v21()
        {
            return (E2 / E1) * v12;
        }
        static double q11()
        {
            return E1 / (1 - (v12 * v21()));
        }
        static double q12()
        {
            return (v12 * E2) / (1 - (v12 * v21()));
        }
        static double q22()
        {
            return E2 / (1 - (v12 * v21()));
        }
        static double q66()
        {
            return G12;
        }

        static double[,] q_matrix = new double[3, 3] { { q11(), q12(), 0 }, { q12(), q22(), 0 }, { 0, 0, q66() } };
        public double t = 0;
        public double top = 0;
        public double bottom = 0;

        public struct stress_and_strains
        {
            public Matrix<double> strains;
            public Matrix<double> stresses;
        }
        public stress_and_strains top_stress_strain = new stress_and_strains();
        public stress_and_strains bottom_stress_strain = new stress_and_strains();

        public double hk1()
        {
            return top - bottom;
        }
        public double hk2()
        {
            return Math.Pow(-top, 2) - Math.Pow(-bottom, 2);
        }
        public double hk3()
        {
            return Math.Pow(top, 3) - Math.Pow(bottom, 3);
        }
        public double q11t()
        {
            return (q11() * Math.Pow(Math.Cos(t), 4)) + 2 * (q12() + 2 * q66()) * (Math.Pow(Math.Sin(t), 2)) * (Math.Pow(Math.Cos(t), 2)) + (q22() * Math.Pow(Math.Sin(t), 4));
        }
        public double q12t()
        {
            return (q11() + q22() - (4 * q66())) * (Math.Pow(Math.Sin(t), 2)) * (Math.Pow(Math.Cos(t), 2)) + (q12() * (Math.Pow(Math.Sin(t), 4) + Math.Pow(Math.Cos(t), 4)));
        }
        public double q22t()
        {
            return (q11() * Math.Pow(Math.Sin(t), 4)) + (2 * (q12() + 2 * q66()) * Math.Pow(Math.Sin(t), 2) * Math.Pow(Math.Cos(t), 2)) + (q22() * Math.Pow(Math.Cos(t), 4));
        }
        public double q16t()
        {
            return ((q11() - q12() - (2 * q66())) * (Math.Sin(t) * Math.Pow(Math.Cos(t), 3))) + ((q12() - q22() + (2 * q66())) * (Math.Pow(Math.Sin(t), 3) * Math.Cos(t)));
        }
        public double q26t()
        {
            return ((q11() - q12() - (2 * q66())) * (Math.Pow(Math.Sin(t), 3) * Math.Cos(t))) + ((q12() - q22() + (2 * q66())) * (Math.Sin(t) * Math.Pow(Math.Cos(t), 3)));
        }
        public double q66t()
        {
            return ((q11() + q22() - (2 * q12()) - (2 * q66())) * (Math.Pow(Math.Sin(t), 2) * Math.Pow(Math.Cos(t), 2))) + (q66() * (Math.Pow(Math.Sin(t), 4) + Math.Pow(Math.Cos(t), 4)));
        }

    }



}


