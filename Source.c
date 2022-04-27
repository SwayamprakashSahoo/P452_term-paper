#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define max 1000

double rates(int n, int U, int M, int A, int B, double H_s, double H_2, double CH3_s, double T_ns, double substrate_temp);

double kmc(double r, double H_s, double H_2, double CH3_s, double T_ns, double substrate_temp);

int substrate[max] = { 0 };
int counts[max];

double rates(int n, int U, int M, int A, int B, double H_s, double H_2, double CH3_s, double T_ns, double substrate_temp)
{
	double R = 1 / (1 + (0.3 * exp(-3430 / substrate_temp) + 0.1 * exp(-4420 / substrate_temp) * (H_2 / H_s)));
	double rate;
	switch (n)
	{
	case 0: //activation of site
		rate = (3.2 * pow(10, (-12)) * pow(substrate_temp, 0.5) * exp(-3430 / substrate_temp)) * H_s * U;
		break;

	case 1: //deactivation of site
		rate = (9.6 * pow(10, (-13)) * (pow(substrate_temp, 0.5)) * H_s + (3.2 * pow(10, (-13)) * pow(substrate_temp, 0.5) * exp(-7850 / substrate_temp) * H_2)) * A;
		break;

	case 3: //rate of adsorption of CH3
		rate = 20 * CH3_s * 3757 * sqrt(T_ns) / (4 * 1.56 * pow(10, 15));
		break;

	case 4: //rate of migration
		rate = (1 - pow((1 - R), 4)) * 6.13 * pow(10, 15) * exp(-128400 / (R * substrate_temp)) * M;
		break;

	case 5: //rate of beta-scission
		rate = (2.083 * pow(10, 10) * substrate_temp * exp(180000 / (R * substrate_temp))) * B;
		break;
	}
	return rate;
}

double kmc(double r, double H_s, double H_2, double CH3_s, double T_ns, double substrate_temp)
{
	double time;
	int grey = 0, magenta = 0, green = 0, red = 0, M = 0, B = 0;
	//count of surface lattice
	for (int i = 0; i < max; i++)
	{
		if (substrate[i] == 0)
			grey++;
		else if (substrate[i] == 1)
			magenta++;
		else if (substrate[i] == 2)
			green++;
		else if (substrate[i] == 3)
		{
			red++;
			if (substrate[i + 1] == 1 || substrate[i - 1] == 1)
				M++;
		}
		if (abs(counts[i] - counts[i - 1] > 1 || abs(counts[i] - counts[i + 1] > 1)))
		{
			B++;
		}

	}

	double rate_sum = 0;
	//if there are n sites

	for (int i = 0; i < max; i++)
	{
		if (substrate[i] == 0)
		{
			rate_sum += rates(0, grey, 0, 0, 0, H_s, H_2, CH3_s, T_ns, substrate_temp); //activation
		}

		else if (substrate[i] == 1)
		{
			rate_sum += rates(0, grey, 0, magenta, 0, H_s, H_2, CH3_s, T_ns, substrate_temp); //deactivation
			rate_sum += rates(0, grey, (red + magenta), 0, 0, H_s, H_2, CH3_s, T_ns, substrate_temp); //CH3 adsorption
		}

		else if (substrate[i] == 2)
		{
			rate_sum += rates(0, green, 0, 0, 0, H_s, H_2, CH3_s, T_ns, substrate_temp); //activation
			rate_sum += rates(0, grey, 0, red, 0, H_s, H_2, CH3_s, T_ns, substrate_temp); //deactivation
		}

		else if (substrate[i] == 3)
		{
			if (substrate[i - 1] == 1 || substrate[i + 1] == 1)
				rate_sum += rates(0, grey, M, red, 0, H_s, H_2, CH3_s, T_ns, substrate_temp); //migration

			if (abs(counts[i] - counts[i - 1] > 1 || abs(counts[i] - counts[i - 1] > 1)))
				rate_sum += rates(0, grey, M, red, B, H_s, H_2, CH3_s, T_ns, substrate_temp); //beta-scission
		}

	}

	time = -(log(r) / rate_sum);
	double temp1 = 0, temp2 = 0;

	//fixing the step
	for (int i = 0; i < max; i++)
	{
		if (substrate[i] == 0)
		{
			temp1 += rates(0, grey, 0, 0, 0, H_s, H_2, CH3_s, T_ns, substrate_temp); //activation

			if ((temp2 / rate_sum) < r < (temp1 / rate_sum))
			{
				substrate[i] = 1;
				break;
			}
			temp2 = temp1;
		}

		if (substrate[i] == 1)
		{
			temp1 += rates(0, grey, 0, magenta, 0, H_s, H_2, CH3_s, T_ns, substrate_temp); //deactivation
			if ((temp2 / rate_sum) < r < (temp1 / rate_sum))
			{
				substrate[i] = 0;
				break;
			}
			temp2 = temp1;
			temp1 += rates(0, grey, (red + magenta), 0, 0, H_s, H_2, CH3_s, T_ns, substrate_temp); //CH3 adsorption
			if ((temp2 / rate_sum) < r < (temp1 / rate_sum))
			{
				counts[i] += 1;
				substrate[i] = 2;
				break;
			}
			temp2 = temp1;
			//formation of side edge
			if (substrate[i - 1] == 2 || substrate[i + 1] == 2 || substrate[i - 1] == 3 && substrate[i + 1] == 3)
			{
				if (substrate[i - 1] == 2 && counts[i - 1] == counts[i])
				{
					substrate[i] = 1;
					substrate[i - 1] = 1;
				}
				else if (substrate[i + 1] == 2 && counts[i + 1] == counts[i])
				{
					substrate[i] = 1;
					substrate[i + 1] = 1;
				}
				if (substrate[i - 1] == 3 && counts[i - 1] == counts[i])
				{
					substrate[i] = 1;
					substrate[i - 1] = 1;
				}
				if (substrate[i + 1] == 3 && counts[i + 1] == counts[i])
				{
					substrate[i] = 1;
					substrate[i + 1] = 1;
				}
			}
		}

		if (substrate[i] == 2)

		{
			temp1 += rates(0, green, 0, 0, 0, H_s, H_2, CH3_s, T_ns, substrate_temp); //activation
			if ((temp2 / rate_sum) < r < (temp1 / rate_sum))
			{
				substrate[i] = 3;
				break;
			}
			temp2 = temp1;
			temp1 += rates(0, grey, 0, red, 0, H_s, H_2, CH3_s, T_ns, substrate_temp); //deactivation
			if ((temp2 / rate_sum) < r < (temp1 / rate_sum))
			{
				counts[i] -= 1;
				substrate[i] = 1;
				break;
			}
			temp2 = temp1;
		}

		if (substrate[i] == 3)
		{
			if (substrate[i - 1] == 1 && substrate[i + 1] == 1)
			{
				if (rand() % 2 == 0)
				{
					temp1 += rates(0, grey, M, red, 0, H_s, H_2, CH3_s, T_ns, substrate_temp); //migrate right
					if ((temp2 / rate_sum) < r < (temp1 / rate_sum))
					{
						counts[i] -= 1;
						counts[i + 1] += 1;
						substrate[i] = 1;
						if (counts[i + 1] < counts[i])
							substrate[i + 1] = 0;
						break;
					}
				}
				else
				{
					temp1 += rates(0, grey, M, red, 0, H_s, H_2, CH3_s, T_ns, substrate_temp); //migrate left
					if ((temp2 / rate_sum) < r < (temp1 / rate_sum))
					{
						counts[i] -= 1;
						counts[i - 1] += 1;
						substrate[i] = 1;
						if (counts[i - 1] < counts[i])
							substrate[i - 1] = 0;
						break;
					}
				}
				temp2 = temp1;
			}
			else if (substrate[i - 1] == 1 && substrate[i + 1] != 1)
			{
				temp1 += rates(0, grey, M, red, 0, H_s, H_2, CH3_s, T_ns, substrate_temp); //migrate left
				if ((temp2 / rate_sum) < r < (temp1 / rate_sum))
				{
					counts[i] -= 1;
					counts[i - 1] += 1;
					substrate[i] = 1;
					if (counts[i - 1] < counts[i])
						substrate[i - 1] = 0;
					break;
				}
				temp2 = temp1;
			}

			else if (substrate[i - 1] != 1 && substrate[i + 1] == 1)
			{
				temp1 += rates(0, grey, M, red, 0, H_s, H_2, CH3_s, T_ns, substrate_temp); //migrate right
				if ((temp2 / rate_sum) < r < (temp1 / rate_sum))
				{
					counts[i] -= 1;
					counts[i + 1] += 1;
					substrate[i] = 1;
					if (counts[i + 1] < counts[i])
						substrate[i + 1] = 0;
					break;
				}
				temp2 = temp1;
			}

			if (abs(counts[i] - counts[i - 1] > 1 || abs(counts[i] - counts[i + 1] > 1)))
			{
				temp1 += rates(0, grey, M, red, B, H_s, H_2, CH3_s, T_ns, substrate_temp); //beta-scission
				if ((temp2 / rate_sum) < r < (temp1 / rate_sum))
				{
					counts[i] -= 2;
					substrate[i] = 1;
					break;
				}
				temp2 = temp1;
			}

			//formation of side edge
			if (substrate[i - 1] == 2 || substrate[i + 1] == 2 || substrate[i - 1] == 3 && substrate[i + 1] == 3)
			{
				if (substrate[i - 1] == 2 && counts[i - 1] == counts[i])
				{
					substrate[i] = 1;
					substrate[i - 1] = 1;
				}
				else if (substrate[i + 1] == 2 && counts[i + 1] == counts[i])
				{
					substrate[i] = 1;
					substrate[i + 1] = 1;
				}
				if (substrate[i - 1] == 3 && counts[i - 1] == counts[i])
				{
					substrate[i] = 1;
					substrate[i - 1] = 1;
				}
				if (substrate[i + 1] == 3 && counts[i + 1] == counts[i])
				{
					substrate[i] = 1;
					substrate[i + 1] = 1;
				}
			}
		}

	}
	return time;
}

//main function to provide parameters like surface temperature, concentration of gases etc.
void main()
{
	double T = 0;
	//creating files for data storage
	FILE* fptr1 = NULL;
	FILE* fptr2 = NULL;
	FILE* fptr3 = NULL;
	FILE* fptr4 = NULL;
	FILE* fptr5 = NULL;
	FILE* fptr6 = NULL;

	int n = 600000;
	//573K
	double H_s = 18 * pow(10, 14), H_2 = 3.37 * pow(10, 17), CH3_s = 2.03 * pow(10, 13), T_ns = 781, substrate_temp = 573;
	fptr1 = fopen("C:/Users/sahoo/Desktop/Computational Asiignments/573K.txt", "w");

	for (int i = 0; i < n; i++)
	{
		T += kmc(rand() / (double)RAND_MAX, H_s, H_2, CH3_s, T_ns, substrate_temp);
	}
	for (int j = 0; j < max; j++)
	{
		fprintf(fptr1, "\t%d\t\n", counts[j]);
	}
	printf("\nSimulation for deposition at 573K complete with time %lf\n", T);
	T = 0;

	//673K
	H_s = 9.34 * pow(10, 14), H_2 = 2.87 * pow(10, 17), CH3_s = 2.50 * pow(10, 13), T_ns = 857, substrate_temp = 673;
	fptr2 = fopen("C:/Users/sahoo/Desktop/Computational Asiignments/673K.txt", "w");

	for (int i = 0; i < n; i++)
	{
		T += kmc(rand() / (double)RAND_MAX, H_s, H_2, CH3_s, T_ns, substrate_temp);
	}
	for (int j = 0; j < max; j++)
	{
		fprintf(fptr2, "\t%d\t\n", counts[j]);
	}
	printf("\nSimulation for deposition at 673K complete with time %lf\n", T);
	T = 0;

	//773K
	H_s = 5.48 * pow(10, 14), H_2 = 2.50 * pow(10, 17), CH3_s = 2.55 * pow(10, 13), T_ns = 935, substrate_temp = 773;
	fptr3 = fopen("C:/Users/sahoo/Desktop/Computational Asiignments/773K.txt", "w");

	for (int i = 0; i < n; i++)
	{
		T += kmc(rand() / (double)RAND_MAX, H_s, H_2, CH3_s, T_ns, substrate_temp);
	}
	for (int j = 0; j < max; j++)
	{
		fprintf(fptr3, "\t%d\t\n", counts[j]);
	}
	printf("\nSimulation for deposition at 773K complete with time %lf\n", T);
	T = 0;

	//973K
	H_s = 2.63 * pow(10, 14), H_2 = 1.99 * pow(10, 17), CH3_s = 2.00 * pow(10, 13), T_ns = 1094, substrate_temp = 973;
	fptr4 = fopen("C:/Users/sahoo/Desktop/Computational Asiignments/973K.txt", "w");

	for (int i = 0; i < n; i++)
	{
		T += kmc(rand() / (double)RAND_MAX, H_s, H_2, CH3_s, T_ns, substrate_temp);
	}
	for (int j = 0; j < max; j++)
	{
		fprintf(fptr4, "\t%d\t\n", counts[j]);
	}
	printf("\nSimulation for deposition at 973K complete with time %lf\n", T);
	T = 0;

	//1173K
	H_s = 1.84 * pow(10, 14), H_2 = 1.65 * pow(10, 17), CH3_s = 1.44 * pow(10, 13), T_ns = 1266, substrate_temp = 1173;
	fptr5 = fopen("C:/Users/sahoo/Desktop/Computational Asiignments/1173K.txt", "w");

	for (int i = 0; i < n; i++)
	{
		T += kmc(rand() / (double)RAND_MAX, H_s, H_2, CH3_s, T_ns, substrate_temp);
	}
	for (int j = 0; j < max; j++)
	{
		fprintf(fptr5, "\t%d\t\n", counts[j]);
	}
	T = 0;

	//1373K
	H_s = 1.63 * pow(10, 14), H_2 = 1.41 * pow(10, 17), CH3_s = 1.13 * pow(10, 13), T_ns = 1442, substrate_temp = 1373;
	fptr6 = fopen("C:/Users/sahoo/Desktop/Computational Asiignments/1373K.txt", "w");

	for (int i = 0; i < n; i++)
	{
		T += kmc(rand() / (double)RAND_MAX, H_s, H_2, CH3_s, T_ns, substrate_temp);
	}
	for (int j = 0; j < max; j++)
	{
		fprintf(fptr6, "\t%d\t\n", counts[j]);
	}
	printf("\nSimulation for deposition at 1373K complete with time %lf\n", T);

	if (fptr1 == NULL || fptr2 == NULL || fptr3 == NULL || fptr4 == NULL || fptr5 == NULL || fptr6 == NULL)
	{
		// File not created hence exit
		printf("\nUnable to create one of the files.\n");
		exit(EXIT_FAILURE);
	}

}
