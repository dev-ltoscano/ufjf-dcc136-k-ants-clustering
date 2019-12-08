#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <limits>
#include <chrono>

using namespace std;
using namespace std::chrono;

// Fator que ajusta a probabilidade de pegar um item
float pickItemFactor = 0.15f;

// Fator que ajusta a probabilidade de soltar um item
float dropItemFactor = 0.15f;

// Fator de dependência entre os dados
float dependencyDataFactor = 0.1f;

/// <summary>
/// Representa um dado para clusterização
/// </summary>
struct DataItem
{
	// Identificador do item
	int id;

	// Tamanho do vetor de dados
	int dataLength;

	// Vetor de dados
	float *data;

	// Distância até o centroide do cluster
	float distToCentroid;

	// Menor distância  até um centroide de um cluster externo
	float distExtCentroid;

	// Indíce de silhueta para o item
	float silhouette;

	// Identificador do cluster associado
	int clusterId;

	DataItem()
	{
		this->id = -1;
		this->dataLength = -1;
		this->data = NULL;

		this->distToCentroid = numeric_limits<float>::max();
		this->distExtCentroid = numeric_limits<float>::max();
		this->silhouette = 0.0f;
	}

	DataItem(int id, int dataLength)
	{
		this->id = id;
		this->dataLength = dataLength;
		this->data = new float[dataLength];

		this->distToCentroid = numeric_limits<float>::max();
		this->distExtCentroid = numeric_limits<float>::max();
		this->silhouette = 0.0f;
	}
};

/// <summary>
/// Representa uma célula da grid bidimensional
/// </summary>
struct GridCell
{
	// Indica se a célula está bloqueada ou não para inserções de itens
	bool isLocked = false;

	// Indica se a célula contém um item
	bool containsItem = false;

	// Item contido na célula
	DataItem item;
};

/// <summary>
/// Representa a posição da célula no grid
/// </summary>
struct GridCellPosition
{
	// Posição (X, Y) da célula no grid
	int cellPosX, cellPosY;

	GridCellPosition()
	{
		this->cellPosX = -1;
		this->cellPosY = -1;
	}

	GridCellPosition(int cellPosX, int cellPosY)
	{
		this->cellPosX = cellPosX;
		this->cellPosY = cellPosY;
	}
};

/// <summary>
/// Representa um grid bidimensional
/// </summary>
struct Grid
{
	// Dimensões do grid bidimensional
	int xDim, yDim;

	// Matriz de células do grid
	GridCell **gridCell;

	// Quantidade de itens armazenados no grid
	int itemCount;

	Grid(int xDim, int yDim)
	{
		this->xDim = xDim;
		this->yDim = yDim;

		this->gridCell = new GridCell*[xDim];

		for (int i = 0; i < xDim; i++)
		{
			this->gridCell[i] = new GridCell[yDim];
		}

		this->itemCount = 0;
	}

	~Grid()
	{
		for (int i = 0; i < xDim; i++)
		{
			delete[] this->gridCell[i];
		}

		delete[] this->gridCell;
	}

	void insert(GridCellPosition position, DataItem item)
	{
		this->gridCell[position.cellPosX][position.cellPosY].item = item;
		this->gridCell[position.cellPosX][position.cellPosY].containsItem = true;
		this->itemCount++;
	}

	void remove(GridCellPosition position)
	{
		this->gridCell[position.cellPosX][position.cellPosY].containsItem = false;
		this->itemCount--;
	}
};

/// <summary>
/// Representa uma formiga para clusterização
/// </summary>
struct Ant
{
	// Posição da formiga no grid
	GridCellPosition position;

	// Lista de itens carregados pela formiga
	unordered_map<int, DataItem> itemList;

	Ant(GridCellPosition position)
	{
		this->position = position;
	}

	/// <summary>
	///	Insere um item na lista de itens da formiga
	/// </summary>
	/// <param name="item">Item a ser inserido</param>
	void insertItem(DataItem item)
	{
		this->itemList.insert(pair<int, DataItem>(item.id, item));
	}

	/// <summary>
	/// Remove um item da lista de itens da formiga
	/// </summary>
	/// <param name="itemId">Identificador do item a ser removido</param>
	void removeItem(int itemId)
	{
		this->itemList.erase(itemId);
	}

	/// <summary>
	/// Retorna a quantidade de elementos carregados pela formiga
	/// </summary>
	/// <returns>Quantidade de elementos carregados pela formiga</returns>
	int itemCount()
	{
		return this->itemList.size();
	}
};

/// <summary>
/// Representa um cluster
/// </summary>
struct Cluster
{
	// Identificador do cluster
	int clusterId = -1;

	// Lista de formigas do cluster
	vector<Ant> antList;

	// Quantidade de itens no cluster
	int itemCount = 0;
};

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
/// Gera um número ponto flutuante aleatório
/// </summary>
/// <returns>Número ponto flutuante aleatório</returns>
float generateRandomFloat()
{
	return ((float)rand() / RAND_MAX);
}

/// <summary>
/// Calcula a distância euclidiana entre dois dados
/// </summary>
/// <param name="first_data_item">Primeiro dado</param>
/// <param name="second_data_item">Segundo dado</param>
/// <param name="dataLength">Tamanho dos dados</param>
/// <returns>Distância euclidiana entre dois dado</returns>
float calculateEuclideanDistance(float *first_data_item, float *second_data_item, int dataLength)
{
	float sumValue = 0.0f;

	for (int i = 0; i < dataLength; i++)
	{
		sumValue += pow((first_data_item[i] - second_data_item[i]), 2);
	}

	return sqrt(sumValue);
}

/// <summary>
/// Calcula o índice no vetor que representa uma matriz triangular inferior de distâncias
/// </summary>
/// <param name="i">Índice i da matriz</param>
/// <param name="j">Índice j da matriz</param>
/// <param name="dimension">Ordem da matriz</param>
/// <returns>Índice no vetor que representa a matriz triangular inferior de distâncias</returns>
int calculateDistanceMatrixIndex(int i, int j, int dimension)
{
	if ((i >= 0) && (i < dimension) && (j >= 0) && (j < dimension))
	{
		if (i == j)
		{
			return 0;
		}
		else if (j > i)
		{
			return ((j * (j - 1) / 2 + i) + 1);
		}
		else
		{
			return ((i * (i - 1) / 2 + j) + 1);
		}
	}

	return -1;
}

/// <summary>
/// Cria um vetor que representa uma matriz triangular inferior de distâncias
/// </summary>
/// <param name="data_vect">Vetor de dados</param>
/// <param name="dataCount">Quantidade de dados</param>
/// <param name="dataLength">Tamanho dos dados</param>
/// <param name="distanceRange"></param>
/// <returns></returns>
float* calculateDistanceMatrix(DataItem *data_vect, int dataCount, int dataLength)
{
	// Quantidade de itens na matriz triangular inferior + 1
	unsigned int triangularLowerMatrixSize = (((dataCount * (dataCount - 1)) / 2) + 1);

	// Vetor que representa a matriz trangular inferior de distâncias
	float *distanceMatrix = new float[triangularLowerMatrixSize];

	// A primeira posição guardará a distância zero
	distanceMatrix[0] = 0.0f;

	float maxDistance = 0.0f;
	float tmpDistance;

	// Calcula a distância entre todos os itens e encontra a distância máxima
	for (int i = 1; i < dataCount; i++)
	{
		for (int j = 0; j < i; j++)
		{
			// Calcula a distância entre o item i e o item j
			tmpDistance = calculateEuclideanDistance(data_vect[i].data, data_vect[j].data, dataLength);

			// Atualiza a distância máxima
			if (tmpDistance > maxDistance)
			{
				maxDistance = tmpDistance;
			}

			distanceMatrix[calculateDistanceMatrixIndex(i, j, dataCount)] = tmpDistance;
		}
	}

	// Normaliza as distâncias
	for (int i = 0; i < triangularLowerMatrixSize; i++)
	{
		distanceMatrix[i] = (distanceMatrix[i] / maxDistance);
	}

	return distanceMatrix;
}

/// <summary>
/// Calcula a função de densidade de vizinhança para um dado 
/// </summary>
/// <param name="cluster_vect">Um cluster</param>
/// <param name="antPerCluster">Quantidade de formigas no cluster</param>
/// <param name="data_item">Dado a ser testado</param>
/// <param name="distanceMatrix">Vetor que representa a matriz de distâncias</param>
/// <param name="distMatrixDimension">Ordem da matriz de distâncias</param>
/// <returns>Valor da função de densidade de vizinhança para o dado</returns>
float calculateNeighborhoodDensityFunction(Cluster *cluster_vect, int antPerCluster, DataItem *data_item, float *distanceMatrix, int distMatrixDimension)
{
	// Valor da função de densidade
	float densityFunctionValue = 0.0f;
	
	// Variáveis auxiliares
	float dissimilarity;
	Ant *tmpAnt;

	for (int k = 0; k < antPerCluster; k++)
	{
		tmpAnt = &cluster_vect->antList[k];

		// Para todos os itens da formiga
		for (unordered_map<int, DataItem>::iterator it = tmpAnt->itemList.begin(); it != tmpAnt->itemList.end(); it++)
		{
			// Obtém a dissimilaridade entre um item a ser inserido/removido do cluster e um item do cluster
			dissimilarity = distanceMatrix[calculateDistanceMatrixIndex(data_item->id, it->first, distMatrixDimension)];

			// Calcula o valor da função de densidade
			densityFunctionValue += (1.0f - (dissimilarity / dependencyDataFactor));
		}

		densityFunctionValue *= (1.0f / 9.0f);
	}

	return max<float>(0.0f, densityFunctionValue);
}

/// <summary>
/// Calcula a probabilidade de um item ser inserido no cluster
/// </summary>
/// <param name="densityFunctionValue">Valor da função de densidade para o item</param>
/// <returns>Probabilidade de inserção do item</returns>
float calculateItemPickProbability(float densityFunctionValue)
{
	if (densityFunctionValue >= pickItemFactor)
	{
		return 1.0f;
	}
	else
	{
		return pow((densityFunctionValue / (pickItemFactor + densityFunctionValue)), 2);
	}
}

/// <summary>
/// Calcula a probabilidade de um item ser removido de um cluster
/// </summary>
/// <param name="densityFunctionValue">Valor da função de densidade para o item</param>
/// <returns>Probabilidade de remoção do item</returns>
float calculateItemDropProbability(float densityFunctionValue)
{
	if (densityFunctionValue >= dropItemFactor)
	{
		return 0.0f;
	}
	else
	{
		return pow((dropItemFactor / (dropItemFactor + densityFunctionValue)), 2);
	}
}

/// <summary>
/// Calcula o índice de silhueta para a clusterização
/// </summary>
/// <param name="cluster_vect">Vetor de clusters</param>
/// <param name="clusterCount">Quantidade de clusters</param>
/// <param name="dataLength">Tamanho dos dados</param>
/// <returns></returns>
float calculateClusterQuality(Cluster *cluster_vect, int clusterCount, int dataLength, DataItem *&centroide_vect)
{
	// Vetor de centroides
	centroide_vect = new DataItem[clusterCount];

	float solutionSilhouette = 0.0f;
	float distCentroid;
	int totalItemCount = 0;

	// Calculando centroides para cada cluster
	for (int k = 0; k < clusterCount; k++)
	{
		// Centroide do cluster k
		centroide_vect[k].id = k;
		centroide_vect[k].dataLength = dataLength;
		centroide_vect[k].data = new float[dataLength];

		// Inicializa os atributos do centroide
		fill_n(centroide_vect[k].data, dataLength, 0.0f);

		// Para cada formiga do cluster k
		for (int i = 0; i < cluster_vect[k].antList.size(); i++)
		{
			// Para cada item da formiga do cluster k
			for (unordered_map<int, DataItem>::iterator it = cluster_vect[k].antList[i].itemList.begin();
				it != cluster_vect[k].antList[i].itemList.end(); it++)
			{
				// Para cada atributo do item da formiga k
				for (int w = 0; w < dataLength; w++)
				{
					// Soma o valor do atributo
					centroide_vect[k].data[w] += it->second.data[w];
				}
			}
		}

		// Para cada atributo do centroide
		for (int j = 0; j < dataLength; j++)
		{
			// Calcula a média de cada atributo do centroide
			centroide_vect[k].data[j] = (centroide_vect[k].data[j] / cluster_vect[k].itemCount);
		}

		// Para cada formiga do cluster k
		for (int i = 0; i < cluster_vect[k].antList.size(); i++)
		{
			// Para cada item da formiga do cluster k
			for (unordered_map<int, DataItem>::iterator it = cluster_vect[k].antList[i].itemList.begin();
				it != cluster_vect[k].antList[i].itemList.end(); it++)
			{
				// Calcula a distância do item até o centroide k
				it->second.distToCentroid = calculateEuclideanDistance(it->second.data, centroide_vect[k].data, dataLength);

				// Marca a distância até um centróide externo como 'infinita'
				it->second.distExtCentroid = numeric_limits<float>::max();
			}
		}
	}

	// Identificador do cluster ao que o item está associado
	int clusterId, clusterExtId;

	// Para cada centroide
	for (int k = 0; k < clusterCount; k++)
	{
		// Id do cluster atual
		int clusterId = k;

		// Para cada formiga do cluster k
		for (int i = 0; i < cluster_vect[k].antList.size(); i++)
		{
			// Para cada item da formiga do cluster k
			for (unordered_map<int, DataItem>::iterator it = cluster_vect[k].antList[i].itemList.begin();
				it != cluster_vect[k].antList[i].itemList.end(); it++)
			{
				// Para todos os centroides
				for (int w = 0; w < clusterCount; w++)
				{
					// Verifica se não é o centroide do cluster atual
					if (w != k)
					{
						// Calcula a distância do item i até o centroide j
						distCentroid = calculateEuclideanDistance(it->second.data, centroide_vect[w].data, dataLength);

						// Verifica se a distância é menor que a atual
						if (distCentroid < it->second.distExtCentroid)
						{
							// Id do cluster mais próximo
							clusterExtId = w;

							// Atualiza a menor distância
							it->second.distExtCentroid = distCentroid;
						}
					}
				}

				// Calcula a silhueta do item i do cluster k
				it->second.silhouette = ((it->second.distExtCentroid - it->second.distToCentroid) / (max(it->second.distExtCentroid, it->second.distToCentroid)));

				// Verifica se o item está no cluster errado
				if (it->second.silhouette < 0)
				{
					// Marca o cluster que o item deveria estar
					it->second.clusterId = clusterExtId;
				}
				else
				{
					// Marca que o item está no cluster correto
					it->second.clusterId = clusterId;
				}

				// Soma a silhueta do item i ao valor de silhueta da solução
				solutionSilhouette += it->second.silhouette;
				totalItemCount++;
			}
		}
	}

	// Retorna a média das silhuetas
	return (solutionSilhouette / totalItemCount);
}

/// <summary>
/// Calcula a quantidade de passos que será feito pela formiga
/// </summary>
/// <param name="currPos">Posição atual da formiga</param>
/// <param name="xDim">Dimensão do X grid</param>
/// <param name="yDim">Dimensão do Y grid</param>
/// <param name="maxStep">Quantidade máxima de passos que a formiga pode dar</param>
/// <returns>Movimento que será feito pela formiga</returns>
GridCellPosition calculateAntMovement(GridCellPosition currPosition, int xDim, int yDim, int maxStep)
{
	// Quantidade mínima de passos
	if (maxStep < 2)
	{
		maxStep = 2;
	}

	// Quantidade aleatória de passos para a formiga
	int stepX = (rand() % maxStep);
	int stepY = (rand() % maxStep);

	// Movimento que será feito pela formiga
	int movementX = (currPosition.cellPosX + stepX);
	int movementY = (currPosition.cellPosY + stepY);

	// Verifica se o movimento ultrapassa os limites do grid
	if (movementX >= xDim)
	{
		movementX = (movementX % xDim);
	}

	if (movementY >= yDim)
	{
		movementY = (movementY % yDim);
	}

	return GridCellPosition(movementX, movementY);
}

/// <summary>
/// Movimenta a formiga pelo grid
/// </summary>
/// <param name="ant">Formiga que será movimentada</param>
/// <param name="grid">Grade bidimensional</param>
/// <param name="maxStep">Máximo de passos que a formiga pode dar</param>
void moveAnt(Ant *ant, Grid *grid, int maxStep)
{
	// Muda a posição da formiga no grid
	ant->position = calculateAntMovement(ant->position, grid->xDim, grid->yDim, maxStep);
}

DataItem* loadDataBin(string fileName, int dataLength)
{
	// Abre o arquivo
	FILE *tensorFile;
	fopen_s(&tensorFile, fileName.c_str(), "rb");

	// Verifica se o arquivo foi aberto
	if (!tensorFile)
	{
		std::cout << "Não foi possível carregar o arquivo " << fileName.c_str() << endl;
		return NULL;
	}

	// Faz a leitura do número de tensores
	int n_tensor;
	fread(&n_tensor, 1, sizeof(int), tensorFile);

	// Vetor de dados
	DataItem *data_vect = new DataItem[n_tensor];

	// Inicializa o vetor de dados
	for (int i = 0; i < n_tensor; i++)
	{
		data_vect[i].id = i;
		data_vect[i].dataLength = dataLength;
		data_vect[i].data = new float[dataLength];
	}

	// Variável auxiliar
	double tensorValue;

	// Faz a leitura dos tensores
	for (int i = 0; i < n_tensor; i++)
	{
		for (int j = 0; j < dataLength; j++)
		{
			fread(&tensorValue, 1, sizeof(double), tensorFile);
			data_vect[i].data[j] = tensorValue;
		}
	}

	// Fecha o arquivo
	fclose(tensorFile);
	return data_vect;
}

DataItem* loadDataText(string fileName, int dataCount, int dataLength)
{
	// Abre o arquivo
	ifstream dataFile(fileName.c_str());

	// Verifica se o arquivo foi aberto
	if (!dataFile.is_open())
	{
		std::cout << "Não foi possível carregar o arquivo " << fileName.c_str() << endl;
		return NULL;
	}

	// Vetor de dados
	DataItem *data_vect = new DataItem[dataCount];

	// Inicializa o vetor de dados
	for (int i = 0; i < dataCount; i++)
	{
		data_vect[i].id = i;
		data_vect[i].dataLength = dataLength;
		data_vect[i].data = new float[dataLength];
	}

	// Faz a leitura dos atributos do item
	for (int i = 0; i < dataCount; i++)
	{
		for (int j = 0; j < dataLength; j++)
		{
			dataFile >> data_vect[i].data[j];
		}
	}

	// Fecha o arquivo
	dataFile.close();
	return data_vect;
}

/// <summary>
/// Salva um arquivo de depuração do algoritmo
/// </summary>
/// <param name="saveFolder">Pasta onde o arquivo será salvo</param>
/// <param name="cluster_vect">Vetor de clusters</param>
/// <param name="clusterCount">Quantidade de clusters</param>
/// <param name="totalIteration">Total de iterações do algoritmo</param>
/// <param name="totalExecutionTime">Tempo total de execução do algoritmo</param>
/// <param name="gridItemCount">Quantidade de itens que sobraram no grid</param>
/// <param name="clusterQuality">Qualidade da clusterização</param>
void saveDebugFile(string saveFolder, Cluster *cluster_vect, int clusterCount, int totalIteration, long totalExecutionTime, int gridItemCount, float clusterQuality)
{
	string fileName = saveFolder + "ant_debug.txt";
	ofstream outFile(fileName.c_str());

	int collectedItemCount = 0;

	for (int i = 0; i < clusterCount; i++)
	{
		outFile << "-- Cluster[" << i << "]: " << endl;

		for (int j = 0; j < cluster_vect[i].antList.size(); j++)
		{
			for (unordered_map<int, DataItem>::iterator it = cluster_vect[i].antList[j].itemList.begin();
				it != cluster_vect[i].antList[j].itemList.end(); it++)
			{
				outFile << "-------- Item[" << it->first << "]" << endl;

			}
		}

		outFile << "-- Cluster item: " << cluster_vect[i].itemCount << endl << endl;
		collectedItemCount += cluster_vect[i].itemCount;
	}

	outFile << "Total iterations: " << totalIteration << endl;
	outFile << "Runtime (Seconds): " << totalExecutionTime << endl;
	outFile << "Total collected items: " << collectedItemCount << endl;
	outFile << "Total remaining items: " << gridItemCount << endl;
	outFile << "Quality of clustering: " << clusterQuality;

	outFile.close();
}

/// <summary>
/// Salva os arquivos do resultado da clusterização
/// </summary>
/// <param name="saveFolder">Pasta onde os arquivos serão salvos</param>
/// <param name="cluster_vect">Vetor de clusters</param>
/// <param name="clusterCount">Quantidade de clusters</param>
/// <param name="cluster_centroid">Vetor dos centroides dos clusters</param>
void saveResultFile(string saveFolder, Cluster *cluster_vect, int clusterCount, DataItem *cluster_centroid, int dataDim)
{
	// Identificador do arquivo do cluster
	int clusterFileId = 1;
	double dataValue;

	// Escreve um arquivo para cada cluster
	for (int i = 0; i < clusterCount; i++)
	{
		// Cria o arquivo do cluster
		string fileName = saveFolder + "cluster" + to_string(clusterFileId) + ".bin";
		FILE *clusterFile;
		fopen_s(&clusterFile, fileName.c_str(), "wb");

		// Verifica se o cluster não foi criado
		if (!clusterFile)
		{
			// Encerra o processo de salvar os arquivos
			std::cout << "Could not create file " << fileName.c_str() << endl;
			break;
		}

		// Escreve a dimensão dos dados
		fwrite(&dataDim, sizeof(int), 1, clusterFile);

		// Escreve todos os dados do centroide do cluster i
		for (int j = 0; j < cluster_centroid[i].dataLength; j++)
		{
			dataValue = static_cast<double>(cluster_centroid[i].data[j]);
			fwrite(&dataValue, sizeof(double), 1, clusterFile);
		}

		// Para cada formiga do cluster i
		for (int j = 0; j < cluster_vect[i].antList.size(); j++)
		{
			// Para cada item da formiga j do cluster i
			for (unordered_map<int, DataItem>::iterator it = cluster_vect[i].antList[j].itemList.begin();
				it != cluster_vect[i].antList[j].itemList.end(); it++)
			{
				// Escreeve todos os dados do item k da formiga j do cluster i
				for (int k = 0; k < it->second.dataLength; k++)
				{
					dataValue = static_cast<double>(it->second.data[k]);
					fwrite(&dataValue, sizeof(double), 1, clusterFile);
				}
			}
		}

		// Fecha o arquivo do cluster
		fclose(clusterFile);

		clusterFileId++;
	}
}

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
	// Verifica se a quantidade mínima de parâmetro foi informada
	if (argc < 9)
	{
		std::cout << "Parâmetros: -[b | t] <filePath> <dataCount> <dataLength> <clusterCount> <antPerCluster> <maxIteration> <outPath>" << endl;
		return 0;
	}

	bool loadFileBinMode = (strcmp(argv[1], "-b") == 0) ? true : false;
	string saveFolder = argv[8];

	if (argc >= 10)
	{
		dependencyDataFactor = atof(argv[9]);
	}

	if (argc >= 12)
	{
		pickItemFactor = atof(argv[10]);
		dropItemFactor = atof(argv[11]);
	}

	bool saveVTK = false;

	if (argc >= 13)
	{
		if (strcmp(argv[12], "-v") == 0)
		{
			saveVTK = true;
		}
	}

	bool quitAfterIteration = false;

	if (argc == 14)
	{
		if (strcmp(argv[13], "-q") == 0)
		{
			quitAfterIteration = true;
		}
	}

	// Define a semente do gerador aleatório
	srand(1);

	// Marca o tempo inicial para calcular o tempo de execução do código
	high_resolution_clock::time_point startTime = high_resolution_clock::now();

	//==================== FASE DE INICIALIZAÇÃO ====================//

	// Quantidade de dados
	int dataCount = atoi(argv[3]);

	// Dimensão e total dos dados
	int dataLength = atoi(argv[4]);

	// Vetor de dados
	std::cout << "Carregando arquivo..." << endl;
	DataItem *data_vect = loadFileBinMode ? loadDataBin(argv[2], dataLength) : loadDataText(argv[2], dataCount, dataLength);

	// Verifica se o arquivo não foi carregado
	if (data_vect == NULL)
	{
		// Encerra o programa
		return 0;
	}

	// Vetor que representa a matriz de distâncias
	std::cout << "Calculando matriz de distâncias..." << endl;
	float *distanceMatrix = calculateDistanceMatrix(data_vect, dataCount, dataLength);

	// Fator que indica o acréscimo de células vazias no grid que serão usadas
	float cellCountFactor = 0.15f;

	// Dimensão do grid bidimensional
	int gridDim = (int)sqrt(((float)dataCount + (cellCountFactor * dataCount)));

	// Grid bidimensional
	Grid grid(gridDim, gridDim);

	// Quantidade máxima inicial de passos que as formigas podem dar
	int maxAntStep = (gridDim / 3);

	// Variáveis auxiliares
	int randomPosX, randomPosY, randomIndex;
	unordered_map<int, GridCellPosition> itemPositionList;
	unordered_map<int, GridCellPosition>::iterator itemPositionListIterator;

	std::cout << "Espalhando os dados na grid..." << endl;

	// Espalha os itens aleatoriamente pelo grid
	for (int i = 0; i < dataCount; i++)
	{
		// Gera uma posição aleatória do grid
		do
		{
			randomPosX = (rand() % gridDim);
			randomPosY = (rand() % gridDim);
		}
		while (grid.gridCell[randomPosX][randomPosY].containsItem);

		GridCellPosition gridCellPosition(randomPosX, randomPosY);

		// Armazena o item na célula do grid
		grid.insert(gridCellPosition, data_vect[i]);

		// Guarda a posição do item no grid
		itemPositionList.insert(pair<int, GridCellPosition>(i, gridCellPosition));
	}

	std::cout << "Quantidade de itens da grid: " << grid.itemCount << endl;

	// Libera a memória do vetor dos dados
	delete[] data_vect;

	// Quantidade de clusters
	int clusterCount = atoi(argv[5]);

	// Quantidade de formigas por cluster
	int antPerCluster = atoi(argv[6]);

	// Vetor de clusters
	Cluster *cluster_vect = new Cluster[clusterCount];

	std::cout << "Espalhando as formigas no grid..." << endl;

	// Espalha as formigas aleatoriamente pelo grid
	for (int i = 0; i < clusterCount; i++)
	{
		// Define o Id do cluster
		cluster_vect[i].clusterId = i;

		// Inicializa a lista de formigas do cluster
		for (int j = 0; j < antPerCluster; j++)
		{
			// Escolhe uma posição aleatória que contenha itens
			randomIndex = (rand() % itemPositionList.size());
			itemPositionListIterator = itemPositionList.begin();
			std::advance(itemPositionListIterator, randomIndex);
			GridCellPosition itemPosition = itemPositionListIterator->second;

			// Remove a posição da lista de posições que contém item
			itemPositionList.erase(itemPositionListIterator);

			// Cria uma formiga e a posiciona no mesmo lugar dos itens
			Ant ant(itemPosition);

			// Insere o item na lista de itens da formiga
			ant.insertItem(grid.gridCell[itemPosition.cellPosX][itemPosition.cellPosY].item);

			// Incrementa a quantidade de itens no cluster
			cluster_vect[i].itemCount++;

			// Remove o item do grid
			grid.remove(itemPosition);

			// Adiciona a formiga ao cluster
			cluster_vect[i].antList.push_back(ant);
		}
	}

	//==================== END FASE DE INICIALIZAÇÃO ====================//


	//==================== FASE DE CLUSTERIZAÇÃO ====================//

	// Número máximo de iterações
	int maxIteration = atoi(argv[7]);

	// Quantidade atual de iterações
	int currIteration = 0;

	// Quantidade de iterações para começar a bloquear células
	int lockIteration = (int)(0.5f * maxIteration);

	// Quantidade de iterações para ocorrer uma redução na velocidade das formigas
	int reduceSpeedIteration = (int)(0.04 * maxIteration);

	// Variáveis auxiliares
	int randomAntIndex;
	float densityFunctionValue, itemProbability, randomFloatNumber;
	Ant *randomAntA, *randomAntB;
	GridCellPosition antPosition;

	std::cout << "Clusterizando..." << endl;

	// Loop principal
	while (currIteration < maxIteration)
	{
		// Para cada cluster
		for (int i = 0; i < clusterCount; i++)
		{
			// Escolhe aleatoriamente uma formiga
			randomAntIndex = (rand() % cluster_vect[i].antList.size());
			randomAntA = &cluster_vect[i].antList[randomAntIndex];

			// Move a formiga de posição no grid
			moveAnt(randomAntA, &grid, maxAntStep);

			// Posição da formiga no grid
			antPosition = randomAntA->position;

			// Verifica se a posição atual da formiga contém itens
			if (grid.gridCell[antPosition.cellPosX][antPosition.cellPosY].containsItem)
			{
				// Item contido na célula do grid
				DataItem gridItem = grid.gridCell[antPosition.cellPosX][antPosition.cellPosY].item;

				// Cálculo da função de densidade de vizinhança para o item
				densityFunctionValue = calculateNeighborhoodDensityFunction(&cluster_vect[i], antPerCluster, &gridItem, distanceMatrix, dataCount);

				// Probabilidade do item ser inserido
				itemProbability = calculateItemPickProbability(densityFunctionValue);

				// Geração de número aleatório
				randomFloatNumber = generateRandomFloat();

				if (randomFloatNumber < itemProbability)
				{
					// Insere o item na lista de itens da formiga
					randomAntA->insertItem(gridItem);

					// Incrementa a quantidade de itens no cluster i
					cluster_vect[i].itemCount++;

					// Remove o item do grid
					grid.remove(antPosition);

					// Remove a posição da lista de posições de itens
					itemPositionList.erase(gridItem.id);

					// Bloqueia a célula após % de iterações
					if (currIteration > lockIteration)
					{
						grid.gridCell[antPosition.cellPosX][antPosition.cellPosY].isLocked = true;
					}
				}
			}
			// Verifica se a formiga tem itens para depositar e se a célula atual não está bloqueada para depositar itens
			else if ((randomAntA->itemList.size() > 0) && !grid.gridCell[antPosition.cellPosX][antPosition.cellPosY].isLocked)
			{
				// Para todos os itens da formiga
				for (unordered_map<int, DataItem>::iterator it = randomAntA->itemList.begin(); 
					it != randomAntA->itemList.end(); it++)
				{
					// Item da formiga
					DataItem antItem = it->second;

					// Cálculo da função de densidade de vizinhança para o item
					densityFunctionValue = calculateNeighborhoodDensityFunction(&cluster_vect[i], antPerCluster, &antItem, distanceMatrix, dataCount);

					// Probabilidade do item ser removido
					itemProbability = calculateItemDropProbability(densityFunctionValue);

					// Geração de número aleatório
					randomFloatNumber = generateRandomFloat();

					if (randomFloatNumber < itemProbability)
					{
						// Insere o item no grid
						grid.insert(antPosition, antItem);

						// Guarda a posição do item no grid
						itemPositionList.insert(pair<int, GridCellPosition>(antItem.id, antPosition));

						// Remove o item da lista de itens da formiga
						it = randomAntA->itemList.erase(it);

						// Decrementa a quantidade de itens no cluster i
						cluster_vect[i].itemCount--;

						// Encerra o loop
						break;
					}
				}
			}

			// Verifica se a formiga tem itens para trocar com outra formiga
			if (randomAntA->itemList.size() > 0)
			{
				// Para todos os clusters
				for (int k = 0; k < clusterCount; k++)
				{
					// Verifica se é um cluster diferente do atual
					if (k != i)
					{
						// Escolhe uma formiga aleatoriamente do cluster k
						randomAntIndex = (rand() % cluster_vect[k].antList.size());
						randomAntB = &cluster_vect[k].antList[randomAntIndex];

						// Para cada item da formiga A
						for (unordered_map<int, DataItem>::iterator it = randomAntA->itemList.begin(); it != randomAntA->itemList.end();)
						{
							// Item da formiga A
							DataItem antItem = it->second;

							// Cálculo da função de densidade de vizinhança para o item
							densityFunctionValue = calculateNeighborhoodDensityFunction(&cluster_vect[k], antPerCluster, &antItem, distanceMatrix, dataCount);

							// Probabilidade do item ser inserido
							itemProbability = calculateItemPickProbability(densityFunctionValue);

							// Geração de número aleatório
							randomFloatNumber = generateRandomFloat();

							// Verifica se o item será inserido na lista de itens da formiga B
							if (randomFloatNumber < itemProbability)
							{
								// Insere o item na lista de itens da formiga B do cluster k
								randomAntB->insertItem(antItem);

								// Incrementa a quantidade de itens no cluster i
								cluster_vect[k].itemCount++;

								// Remove o item na lista de itens da formiga A do cluster i
								it = randomAntA->itemList.erase(it);

								// Decrementa a quantidade de itens no cluster i
								cluster_vect[i].itemCount--;
							}
							else
							{
								it++;
							}
						}
					}
				}
			}

			// Verifica se a formiga está sem itens e ainda há itens no grid
			if ((randomAntA->itemList.size() == 0) && (itemPositionList.size() > 0))
			{
				// Escolhe uma posição que contenha itens para a formiga
				randomIndex = (rand() % itemPositionList.size());
				itemPositionListIterator = itemPositionList.begin();
				std::advance(itemPositionListIterator, randomIndex);
				GridCellPosition itemPosition = itemPositionListIterator->second;

				// Movimenta a formiga para a posição
				randomAntA->position = itemPosition;

				// Insere o item na lista de itens da formiga
				randomAntA->insertItem(grid.gridCell[itemPosition.cellPosX][itemPosition.cellPosY].item);

				// Incrementa a quantidade de itens no cluster i
				cluster_vect[i].itemCount++;

				// Remove o item do grid
				grid.remove(itemPosition);

				// Remove a posição da lista de posições que contém item
				itemPositionList.erase(itemPositionListIterator);
			}
		}

		// Verifica se a iteração atual é a de redução da velocidade das formigas
		if ((currIteration % reduceSpeedIteration) == 0)
		{
			// Reduz a velocidade máxima das formigas
			maxAntStep--;
		}

		if ((currIteration % 100) == 0)
		{
			std::cout << "Iteration: " << currIteration << " | Remain itens: " << grid.itemCount << endl;
		}

		currIteration++;
	}

	// Verifica se ainda sobrou itens no grid após o número máximo de iterações
	if (!quitAfterIteration && grid.itemCount > 0)
	{
		std::cout << "Alocando os itens que sobraram no grid para o melhor cluster provável..." << endl;

		float maxItemProbability = 0.0f;
		int clusterInsertId = -1;

		// Para cada item que restou no grid
		itemPositionListIterator = itemPositionList.begin();

		while (itemPositionListIterator != itemPositionList.end())
		{
			// Posição do item
			GridCellPosition itemPosition = itemPositionListIterator->second;

			// Item contido no grid
			DataItem gridItem = grid.gridCell[itemPosition.cellPosX][itemPosition.cellPosY].item;

			for (int j = 0; j < clusterCount; j++)
			{
				// Cálculo da função de densidade de vizinhança para o item
				densityFunctionValue = calculateNeighborhoodDensityFunction(&cluster_vect[j], antPerCluster, &gridItem, distanceMatrix, dataCount);

				// Probabilidade do item ser inserido no cluster j
				itemProbability = calculateItemPickProbability(densityFunctionValue);

				if (itemProbability >= maxItemProbability)
				{
					clusterInsertId = j;
					maxItemProbability = itemProbability;
				}
			}

			// Escolhe aleatoriamente uma formiga do cluster que terá o item inserido
			randomAntIndex = (rand() % cluster_vect[clusterInsertId].antList.size());
			randomAntA = &cluster_vect[clusterInsertId].antList[randomAntIndex];

			// Insere o item na lista de itens da formiga do cluster que terá o item inserido
			randomAntA->insertItem(gridItem);

			// Incrementa a quantidade de itens no cluster i
			cluster_vect[clusterInsertId].itemCount++;

			// Remove o item do grid
			grid.remove(itemPosition);

			// Remove o item da lista de posições
			itemPositionListIterator = itemPositionList.erase(itemPositionListIterator);
		}
	}

	// Vetor de centroides dos clusters
	DataItem *centroide_vect = NULL;

	// Faz o cálculo do índice de silhueta da solução atual
	float currClusterQuality = calculateClusterQuality(cluster_vect, clusterCount, dataLength, centroide_vect);
	float finalClusterQuality;

	// Quantidade máxima de iterações da fase de refinamento
	int maxLocalSearchIteration = 1000;
	int currLocalSearchIteration = 0;

	std::cout << "Fazendo refinamento da solução..." << endl;
	
	do
	{
		// Atualiza o índice de silhueta final
		finalClusterQuality = currClusterQuality;

		// Desaloca o vetor de centroides se necessário
		if (centroide_vect != NULL)
		{
			delete[] centroide_vect;
		}

		// Faz refinamento na solução
		for (int i = 0; i < clusterCount; i++)
		{
			for (int j = 0; j < antPerCluster; j++)
			{
				for (unordered_map<int, DataItem>::iterator it = cluster_vect[i].antList[j].itemList.begin();
					it != cluster_vect[i].antList[j].itemList.end();)
				{
					// Verifica se o item está no cluster errado
					if (it->second.clusterId != i)
					{
						// Seleciona uma formiga aleatória do cluster que a formiga deseja estar
						randomAntIndex = (rand() % cluster_vect[it->second.clusterId].antList.size());
						randomAntB = &cluster_vect[it->second.clusterId].antList[randomAntIndex];

						// Insere o item no cluster mais adequado
						randomAntB->insertItem(it->second);

						// Incrementa a quantidade de itens no cluster que recebeu o item
						cluster_vect[it->second.clusterId].itemCount++;

						// Remove o item na lista de itens da formiga j do cluster i
						randomAntA = &cluster_vect[i].antList[j];
						it = randomAntA->itemList.erase(it);

						// Decrementa a quantidade de itens no cluster i
						cluster_vect[i].itemCount--;
					}
					else
					{
						it++;
					}
				}
			}
		}

		// Faz o cálculo do índice de silhueta para a nova solução
		currClusterQuality = calculateClusterQuality(cluster_vect, clusterCount, dataLength, centroide_vect);
		
		currLocalSearchIteration++;
	} 
	while ((currClusterQuality > finalClusterQuality) && (currLocalSearchIteration < maxLocalSearchIteration));

	// Marca o tempo final para calcular o tempo de execução do código
	high_resolution_clock::time_point endTime = high_resolution_clock::now();

	// Calcula o tempo total de execução em segundos
	long totalExecutionTime = duration_cast<seconds>(endTime - startTime).count();

	std::cout << "Salvando arquivos finais..." << endl;

	// Salva arquivo de depuração do algoritmo
	saveDebugFile(saveFolder, cluster_vect, clusterCount, currIteration, totalExecutionTime, grid.itemCount, finalClusterQuality);

	int dataDim = (loadFileBinMode) ? sqrt(dataLength) : dataLength;
	saveResultFile(saveFolder, cluster_vect, clusterCount, centroide_vect, dataDim);

	if (saveVTK)
	{
		if (dataLength > 3)
		{
			std::cout << "Não é possível visualizar esse tipo de dado. [ dataLength > 3 ]" << endl;
		}
		else
		{
			VTK vtk;

			float centroidColor = 1.0f;
			float dataItemColor = 100000.0f;

			for (int i = 0; i < clusterCount; i++)
			{
				// Centroide do cluster
				if (dataLength == 1)
				{
					vtk.pointVect.push_back(Point3D(centroide_vect[i].data[0], 0.0f, 0.0f));
				}
				else if (dataLength == 2)
				{
					vtk.pointVect.push_back(Point3D(centroide_vect[i].data[0], centroide_vect[i].data[1], 0.0f));
				}
				else
				{
					vtk.pointVect.push_back(Point3D(centroide_vect[i].data[0], centroide_vect[i].data[1], centroide_vect[i].data[2]));
				}

				vtk.scalarVect.push_back(centroidColor);

				// Itens do cluster
				for (int j = 0; j < antPerCluster; j++)
				{
					for (unordered_map<int, DataItem>::iterator it = cluster_vect[i].antList[j].itemList.begin();
						it != cluster_vect[i].antList[j].itemList.end(); it++)
					{
						if (dataLength == 1)
						{
							vtk.pointVect.push_back(Point3D(it->second.data[0], 0.0f, 0.0f));
						}
						else if (dataLength == 2)
						{
							vtk.pointVect.push_back(Point3D(it->second.data[0], it->second.data[1], 0.0f));
						}
						else
						{
							vtk.pointVect.push_back(Point3D(it->second.data[0], it->second.data[1], it->second.data[2]));
						}

						vtk.scalarVect.push_back(dataItemColor);
					}
				}

				dataItemColor += 100000.0f;
			}

			saveVTKFile(saveFolder, vtk);
		}
	}

	// Desalocando memória
	delete[] centroide_vect;
	delete[] distanceMatrix;
	delete[] cluster_vect;

	return 0;
}