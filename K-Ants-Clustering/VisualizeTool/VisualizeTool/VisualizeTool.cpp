#include <iostream>
#include <string>
#include <filesystem>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace std::experimental::filesystem;

struct Point3D
{
	float x, y, z;

	Point3D(float x, float y, float z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}
};

struct VTK
{
	vector<Point3D> pointVect;
	vector<float> scalarVect;
};

/// <summary>
/// Salva um arquivo no formato VTK para visualização no programa Paraview
/// </summary>
/// <param name="saveFolder">Pasta onde o arquivo será salvo</param>
/// <param name="vtk">Arquivo VTK</param>
void saveVTKFile(string saveFolder, VTK vtk)
{
	// Cria o arquivo informado
	string fileName = saveFolder + "visualize.vtk";
	ofstream vtkFile(fileName.c_str());

	// Verifica se o arquivo foi aberto
	if (!vtkFile.is_open())
	{
		// Encerra o processo de salvar os arquivos
		std::cout << "Could not create file " << fileName.c_str() << endl;
		return;
	}

	vtkFile << std::fixed << std::setprecision(7);

	// Cabeçalho do arquivo VTK
	vtkFile << "# vtk DataFile Version 3.0" << endl;
	vtkFile << "VTK Output" << endl;
	vtkFile << "ASCII" << endl;
	vtkFile << "DATASET POLYDATA" << endl;
	vtkFile << "POINTS  " << vtk.pointVect.size() << "  float" << endl;

	// Pontos do arquivo VTK
	for (int i = 0; i < vtk.pointVect.size(); i++)
	{
		vtkFile << vtk.pointVect[i].x << "  " << vtk.pointVect[i].y << "  " << vtk.pointVect[i].z << endl;
	}

	vtkFile << "POINT_DATA " << vtk.pointVect.size() << endl;
	vtkFile << "SCALARS color float 1" << endl;
	vtkFile << "LOOKUP_TABLE default" << endl;

	// Scalar do arquivo VTK
	for (int i = 0; i < vtk.scalarVect.size(); i++)
	{
		vtkFile << vtk.scalarVect[i] << endl;
	}

	vtkFile.close();
}

int main(int argc, char **argv)
{
	// Verifica se o parâmetro foi informado
	if (argc < 2)
	{
		std::cout << "Parâmetros: <folderPath>" << endl;
		return 0;
	}

	// Caminho da pasta contendo os arquivos para leitura
	string folderPath = argv[1];
	
	// Dimensão dos dados
	int dataDim;
	double xValue, yValue, zValue;

	// Arquivo VTK
	VTK vtkFile;

	float centroidScalar = 1.0f;
	float dataScalar = 100000.0f;

	for (auto &p : directory_iterator(folderPath))
	{
		FILE *clusterFile = fopen(p.path().string().c_str(), "rb");
		fread(&dataDim, 1, sizeof(int), clusterFile);

		// Verifica se o arquivo foi aberto
		if (!clusterFile)
		{
			std::cout << "Não foi possível carregar o arquivo " << p << endl;
			break;
		}
		
		if (dataDim > 3)
		{
			std::cout << "Não é possível visualizar este tipo de dados. [ dataDim > 3 ]" << endl;
			break;
		}

		// Lendo o centroide do cluster
		if (dataDim == 1)
		{
			fread(&xValue, 1, sizeof(double), clusterFile);

			vtkFile.pointVect.push_back(Point3D(xValue, 0.0f, 0.0f));
			vtkFile.scalarVect.push_back(centroidScalar);

			while (fread(&xValue, 1, sizeof(double), clusterFile))
			{
				vtkFile.pointVect.push_back(Point3D(xValue, 0.0f, 0.0f));
				vtkFile.scalarVect.push_back(dataScalar);
			}
		}
		else if (dataDim == 2)
		{
			fread(&xValue, 1, sizeof(double), clusterFile);
			fread(&yValue, 1, sizeof(double), clusterFile);

			vtkFile.pointVect.push_back(Point3D(xValue, yValue, 0.0f));
			vtkFile.scalarVect.push_back(centroidScalar);

			while (fread(&xValue, 1, sizeof(double), clusterFile) && fread(&yValue, 1, sizeof(double), clusterFile))
			{
				vtkFile.pointVect.push_back(Point3D(xValue, yValue, 0.0f));
				vtkFile.scalarVect.push_back(dataScalar);
			}
		}
		else
		{
			fread(&xValue, 1, sizeof(double), clusterFile);
			fread(&yValue, 1, sizeof(double), clusterFile);
			fread(&zValue, 1, sizeof(double), clusterFile);

			vtkFile.pointVect.push_back(Point3D(xValue, yValue, zValue));
			vtkFile.scalarVect.push_back(centroidScalar);

			while (fread(&xValue, 1, sizeof(double), clusterFile) && fread(&yValue, 1, sizeof(double), clusterFile) && fread(&zValue, 1, sizeof(double), clusterFile))
			{
				vtkFile.pointVect.push_back(Point3D(xValue, yValue, zValue));
				vtkFile.scalarVect.push_back(dataScalar);
			}
		}

		fclose(clusterFile);
		dataScalar += 100000.0f;
	}

	saveVTKFile(folderPath, vtkFile);
	std::cout << "OK!" << endl;

	return 0;
}