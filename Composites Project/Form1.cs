using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Composites_Project
{
    public partial class Form1 : Form
    {
        List<compositeLayer> layers = new List<compositeLayer>();           //Create list to store all the composite layers.
        List<TextBox> angles = new List<TextBox>();                         //Create list to hold all the angles of the plys that are entered into the textboxes.
        int layerToDisplay = 0;                                             //Used to display a certain layer's info on the form.
        Matrix<double> midplane_strains;                                    //Define the midplane strains matrix.
        Matrix<double> midplane_curvatures;                                 //Define the midplane curvatures matrix.


        public Form1()
        {
            /*Load all the form components.*/
            InitializeComponent();

            /*Add all the angle textboxes to the list.*/
            angles.Add(layer1Angle_textBox);
            angles.Add(layer2Angle_textBox);
            angles.Add(layer3Angle_textBox);
            angles.Add(layer4Angle_textBox);
            angles.Add(layer5Angle_textBox);
            angles.Add(layer6Angle_textBox);
        }

        /*Event handler runs anytime of the data textboxes changes.*/
        private void textBoxTextChanged(object sender, EventArgs e)
        {
            /*Start by clearing all the previous layers information.*/
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
                double top = (compositeLayer.numberOfLayers * compositeLayer.layerThickness) / 2;
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
                    /*Do nothing on an error.*/
                    catch { }

                    /* Determine the top and bottom of each layer. */
                    layers[i].top = top - (i * compositeLayer.layerThickness);
                    layers[i].bottom = top - ((i + 1) * compositeLayer.layerThickness);
                }

                compositeLayer.processTemp = Convert.ToDouble(processTemp_textBox.Text);
                compositeLayer.roomTemp = Convert.ToDouble(roomTemp_textBox.Text);
                compositeLayer.alpha_1 = Convert.ToDouble(alpha1_textBox.Text);
                compositeLayer.alpha_2 = Convert.ToDouble(alpha2_textBox.Text);
            }
            /*If any of the textboxes doesn't have a number entered in it skip the rest of the function and do nothing.*/
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
            
            /*deltaT is the change in temperature.*/
            double deltaT = compositeLayer.roomTemp - compositeLayer.processTemp; 
            
            /*Calculate the A,B,D matrix values.*/
            for (int i = 0; i < layers.Count; i++)
            {
                q11a += layers[i].q11t() * layers[i].hk1();
                q11b += layers[i].q11t() * (layers[i].hk2());
                q11d += layers[i].q11t() * (layers[i].hk3());

                q12a += layers[i].q12t() * layers[i].hk1();
                q12b += layers[i].q12t() * (layers[i].hk2());
                q12d += layers[i].q12t() * (layers[i].hk3());

                q16a += layers[i].q16t() * layers[i].hk1();
                q16b += layers[i].q16t() * (layers[i].hk2());
                q16d += layers[i].q16t() * (layers[i].hk3());

                q26a += layers[i].q26t() * layers[i].hk1();
                q26b += layers[i].q26t() * (layers[i].hk2());
                q26d += layers[i].q26t() * (layers[i].hk3());

                q22a += layers[i].q22t() * layers[i].hk1();
                q22b += layers[i].q22t() * (layers[i].hk2());
                q22d += layers[i].q22t() * (layers[i].hk3());

                q66a += layers[i].q66t() * layers[i].hk1();
                q66b += layers[i].q66t() * (layers[i].hk2());
                q66d += layers[i].q66t() * (layers[i].hk3());



                /***** Calculate thermal strains ******/
                /*Create a new matrix for thermal strains.*/
                Matrix<double> thermal_strains = DenseMatrix.OfArray(new double[3, 1]);

                /*Calculate alpha x, y, and xy.*/
                layers[i].alpha_x = (compositeLayer.alpha_1 * Math.Pow(Math.Cos(layers[i].t), 2)) + (compositeLayer.alpha_2 * Math.Pow(Math.Sin(layers[i].t), 2));
                layers[i].alpha_y = (compositeLayer.alpha_1 * Math.Pow(Math.Sin(layers[i].t), 2)) + (compositeLayer.alpha_2 * Math.Pow(Math.Cos(layers[i].t), 2));
                layers[i].alpha_xy = 2 * Math.Cos(layers[i].t) * Math.Sin(layers[i].t) * (compositeLayer.alpha_1 - compositeLayer.alpha_2);
                
                /*Populate the thermal strains matrix.*/
                thermal_strains[0, 0] = layers[i].alpha_x * deltaT;
                thermal_strains[1,0] = layers[i].alpha_y * deltaT;
                thermal_strains[2, 0] = layers[i].alpha_xy * deltaT;

                /*Assign the thermal strains matrix to the compositeLayer object.*/
                layers[i].thermal_strains.strains = thermal_strains;
            }

            /*Assemble the ABD matrix(cies).*/
            Matrix<double> a_matrix = DenseMatrix.OfArray(new double[,]{ {q11a,q12a,q16a, },
                                                                        {q12a,q22a,q26a },
                                                                        {q16a,q26a,q66a } });
            Matrix<double> b_matrix = DenseMatrix.OfArray(new double[,]{ { q11b / 2, q12b / 2, q16b / 2, },
                                                                        { q12b / 2, q22b / 2, q26b / 2 },
                                                                        { q16b / 2, q26b / 2, q66b / 2 } });
            Matrix<double> d_matrix = DenseMatrix.OfArray(new double[,]{ { q11d / 3, q12d / 3, q16d / 3, },
                                                                        { q12d / 3, q22d / 3, q26d / 3 },
                                                                        { q16d / 3, q26d / 3, q66d / 3 } });
            Matrix<double> abd_matrix = DenseMatrix.OfArray(new double[,]{ { q11a, q12a, q16a, q11b / 2, q12b / 2, q16b / 2 },
                                                                            { q12a, q22a, q26a, q12b / 2, q22b / 2, q26b / 2 },
                                                                            { 2 * q16a, 2 * q26a, 2 * q66a, q16b / 2, q26b / 2, q66b / 2 },
                                                                            { q11b / 2, q12b / 2, q16b / 2, q11d / 3, q12d / 3, q16d / 3 },
                                                                            { q12b / 2, q22b / 2, q26b / 2, q12d / 3, q22d / 3, q26d / 3 },
                                                                            { q16b / 2, q26b / 2, q66b / 2, q16d / 3, q26d / 3, q66d / 3 } });
           /*Get the inverse ABD matrix.*/
            Matrix<double> abd_inv_matrix = abd_matrix.Inverse();

            /*Find the thermal force resultant matrix.*/
            double Nxt = 0;
            double Nyt = 0;
            double Nxyt = 0;
            double Mxt = 0;
            double Myt = 0;
            double Mxyt = 0;

            for (int i = 0; i < compositeLayer.numberOfLayers; i++)
            {
                Nxt += ((layers[i].q11t() * layers[i].alpha_x) + (layers[i].q12t() * layers[i].alpha_y) + (layers[i].q16t() * layers[i].alpha_xy)) * layers[i].hk1();
                Nyt += ((layers[i].q12t() * layers[i].alpha_x) + (layers[i].q22t() * layers[i].alpha_y) + (layers[i].q26t() * layers[i].alpha_xy)) * layers[i].hk1();
                Nxyt += ((layers[i].q16t() * layers[i].alpha_x) + (layers[i].q26t() * layers[i].alpha_y) + (layers[i].q66t() * layers[i].alpha_xy)) * layers[i].hk1();
                Mxt += ((layers[i].q11t() * layers[i].alpha_x) + (layers[i].q12t() * layers[i].alpha_y) + (layers[i].q16t() * layers[i].alpha_xy)) * layers[i].hk2();
                Myt += ((layers[i].q12t() * layers[i].alpha_x) + (layers[i].q22t() * layers[i].alpha_y) + (layers[i].q26t() * layers[i].alpha_xy)) * layers[i].hk2();
                Mxyt += ((layers[i].q16t() * layers[i].alpha_x) + (layers[i].q26t() * layers[i].alpha_y) + (layers[i].q66t() * layers[i].alpha_xy)) * layers[i].hk2();
            }
            Nxt *= deltaT;
            Nyt *= deltaT;
            Nxyt *= deltaT;
            Mxt *= deltaT / 2;
            Myt *= deltaT / 2;
            Mxyt *= deltaT / 2;

            /*Assemble the thermal force resultant matrix.*/
            Matrix<double> thermal_force_matrix = DenseMatrix.OfArray(new double[,] { { Nxt },
                                                                                        { Nyt },
                                                                                        { Nxyt },
                                                                                        { Mxt },
                                                                                        { Myt },
                                                                                        { Mxyt } });

            /*Assemble force matrix due to external forces.*/
            Matrix<double> forces_matrix = DenseMatrix.OfArray(new double[,] { { compositeLayer.Nx },
                                                                                { compositeLayer.Ny },
                                                                                { compositeLayer.Nxy },
                                                                                { compositeLayer.Mx },
                                                                                { compositeLayer.My },
                                                                                { compositeLayer.Mxy } });
            
            /*Add the thermal forces to the external forces to complete the force resultant matrix.*/
            forces_matrix = forces_matrix.Add(thermal_force_matrix);

            /*Solve for the strain and curvature matrix of the whole composite part.*/
            Matrix<double> strain_and_curvature = abd_inv_matrix.Multiply(forces_matrix);
            
            /*Break the strain and curvature matrix apart into two seperate matricies, one for strain, one for curvatures.*/
            midplane_strains = DenseMatrix.OfArray(new double[,] { { strain_and_curvature[0,0]},
                                                                    { strain_and_curvature[1,0] },
                                                                    { strain_and_curvature[2,0] }});
            midplane_curvatures = DenseMatrix.OfArray(new double[,] { { strain_and_curvature[3, 0] },
                                                                    { strain_and_curvature[4, 0] },
                                                                    { strain_and_curvature[5, 0] }});

            /*Calculate the stress and strain for each layer.*/
            for (int i = 0; i < compositeLayer.numberOfLayers; i++)
            {
                Matrix<double> qBarMatrix = DenseMatrix.OfArray(new double[,] { { layers[i].q11t(), layers[i].q12t(), layers[i].q16t()},
                                                                                { layers[i].q12t(), layers[i].q22t(), layers[i].q26t()},
                                                                                { layers[i].q16t(), layers[i].q26t(), layers[i].q66t()},});
                layers[i].top_stress_strain.strains = midplane_strains.Add(midplane_curvatures.Multiply(layers[i].top));
                layers[i].bottom_stress_strain.strains = midplane_strains.Add(midplane_curvatures.Multiply(layers[i].bottom));
                layers[i].top_stress_strain.stresses = qBarMatrix.Multiply(layers[i].top_stress_strain.strains.Subtract(layers[i].thermal_strains.strains));
                layers[i].bottom_stress_strain.stresses = qBarMatrix.Multiply(layers[i].bottom_stress_strain.strains.Subtract(layers[i].thermal_strains.strains));
                displayLayer(0);
            }
        }

        /*displayLayer function displays a certain layers stress, strain, and Q bar matrix in the textboxes on the form.*/
        void displayLayer(int layerNumber)
        {
            try
            {
                compositeLayer thisLayer = layers[layerNumber];                
                int roundingDigits = 10;
                q11t_textBox.Text = Convert.ToString(Math.Round(thisLayer.q11t(), roundingDigits));
                q12t_textBox.Text = Convert.ToString(Math.Round(thisLayer.q12t(), roundingDigits));
                q12t_textBox2.Text = Convert.ToString(Math.Round(thisLayer.q12t(), roundingDigits));
                q16t_textBox.Text = Convert.ToString(Math.Round(thisLayer.q16t(), roundingDigits));
                q16t_textBox2.Text = Convert.ToString(Math.Round(thisLayer.q16t(), roundingDigits));
                q22t_textBox.Text = Convert.ToString(Math.Round(thisLayer.q22t(), roundingDigits));
                q26t_textBox.Text = Convert.ToString(Math.Round(thisLayer.q26t(), roundingDigits));
                q26t_textBox2.Text = Convert.ToString(Math.Round(thisLayer.q26t(), roundingDigits));
                q66t_textBox.Text = Convert.ToString(Math.Round(thisLayer.q66t(), roundingDigits));

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
            displayedLayer_label.Text = "Layer " + Convert.ToString(layerToDisplay + 1);
        }

        /*nextLayerButton click increments layeToDisplay when the next layer button is clicked.*/
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

        /*previousLayerButton click decrements layerToDisplay when the previous layer button is clicked.*/
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

    }


    /*compositeLayer is the main class that holds all relevant information to any given layer with the composite.*/
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
        public static double processTemp;
        public static double roomTemp;
        public static double alpha_1;
        public static double alpha_2;

        public double alpha_x;
        public double alpha_y;
        public double alpha_xy;
        public double beta_x;
        public double beta_y;
        public double beta_xy;

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
        public stress_and_strains thermal_strains = new stress_and_strains();
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


