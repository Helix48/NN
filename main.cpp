#include <iostream>
#include <vector>
#include <cstdlib>
#include <cassert>
#include <cmath>

using namespace std;

struct Connection
{
    double weight;
    double deltaWeight;
};

class Neuron;

typedef vector<Neuron> Layer;

// ****************** class Neuron ********************

class Neuron
{
public:
    Neuron(unsigned numOutputs, unsigned myIndex);
    void setOutputVal(double val) { m_outputVal=val; }
    double getOutputVal(void) const { return m_outputVal; }
    void feedForward(const Layer &prevLayer);

private:
    static double transferFunction(double x);
    static double transferFunctionDerivative(double x);
    static double randomWeight(void) { return rand() / double{RAND_MAX}; } //generates a random number from 0-1
    double m_outputVal;
    vector<Connection> m_outputWeights;
    unsigned m_myIndex;
};

double Neuron::transferFunction(double x)
{
    //tanh - output range [-1.0..1.0]
    return tanh(x);
}

double Neuron::transferFunctionDerivative(double x)
{
    // tanh derivative
    return 1.0 - x * x;
}

void Neuron::feedForward(const Layer &prevLayer)
{
    double sum=0.0;

    // Sum the previous layer's outputs (which are our inputs)
    // Include the bias node from the previous layer.

    for (unsigned n=0;n<prevLayer.size();++n)
    {
        sum+=prevLayer[n].getOutputVal() *
            prevLayer[n].m_outputWeights[m_myIndex].weight;
    }

    m_outputVal = Neuron::transferFunction(sum);
}

Neuron::Neuron(unsigned numOutputs, unsigned myIndex)
{
    for (unsigned c=0;c<numOutputs; ++c)
    {
        m_outputWeights.push_back(Connection());
        m_outputWeights.back().weight = randomWeight();
    }

    m_myIndex = myIndex;
}

// ****************** class Net ***********************

class Net
{
public:
    Net(const vector<unsigned> &topology);
    void feedForward(const vector<double> &inputVals);
    void backProp(const vector<double> &targetVals);
    void getResults(vector<double> &resultVals) const {};

private:
    vector<Layer> m_layers; // m_layers[layerNum][neuronNum]
    double m_error;
    double m_recentAverageError;
    double m_recentAverageSmoothingFactor;
};

void Net::backProp(const vector<double> &targetVals)
{
    //Calculate overall net error (RMS of output neuron errors)

    Layer &outputLayer = m_layers.back();
    m_error = 0.0;

    for (unsigned n=0;n<outputLayer.size()-1;++n)
    {
        double delta = targetVals[n] - outputLayer[n].getOutputVal();
    }
    m_error /= outputLayer.size() - 1; //get average error squared
    m_error = sqrt(m_error); // RMS

    /*
    //Implement a recent average measurement:

    m_recentAverageError=
    (m_recentAverageError*m_recentAverageSmoothingFactor+m_error)
    / (m_recentAverageSmoothingFactor+1.0);
    */

    //Calculate output layer gradients

    for (unsigned n=0;n<outputLayer.size()-1;++n)
    {
        outputLayer[n].calcOutputGradients(targetVals[n]);
    }

    //Calculate gradients on hidden layers

    for (unsigned layerNum=m_layers.size()-2;layerNum>0;--layerNum)
    {
        Layer &hiddenLayer = m_layers[layerNum];
        Layer &nextLayer = m_layers[layerNum+1];

        for (unsigned n=0;n<hiddenLayer.size();++n)
        {
            hiddenLayer[n].calcHiddenGradients(nextLayer);
        }
    }

    //For all layers from outputs to first hidden layer,
    //update connection weights
}

void Net::feedForward(const vector<double> &inputVals)
{
    assert(inputVals.size() == m_layers[0].size() - 1);

    //Assign (latch) the input values into the input neurons
    for (unsigned i=0;i<inputVals.size();++i)
    {
        m_layers[0][i].setOutputVal(inputVals[i]);
    }

    //Forward propagate
    for (unsigned layerNum=1;layerNum<m_layers.size();++layerNum)
    {
        Layer &prevLayer = m_layers[layerNum-1];
        for (unsigned n=0;n<m_layers[layerNum].size() - 1;++n)
        {
            m_layers[layerNum][n].feedForward();
        }
    }
}

Net::Net(const vector<unsigned> &topology)
{
    unsigned numLayers=topology.size();
    for (unsigned layerNum=0;layerNum<numLayers; ++layerNum)
    {
        m_layers.push_back(Layer());
        unsigned numOutputs = layerNum == topology.size() - 1 ? 0 : topology[layerNum + 1]; //if the loop is on the output layer (last layer) outputs = 0 if not outputs = amount of nodes to feed into next layer

        //We have made a new layer, now to fill it with neurons, and
        //add a bias neuron to the layer:
        for (unsigned neuronNum=0;neuronNum<=topology[layerNum];++neuronNum)
        {
            m_layers.back().push_back(Neuron(numOutputs, neuroNum));
            cout<<"Made a Neuron!"<<endl;
        }
    }
}

int main()
{
    // e.g., { 3, 2, 1 }
    vector<unsigned> topology;
    topology.push_back(3);
    topology.push_back(2);
    topology.push_back(1);
    Net myNet(topology);

    vector<double> inputVals;
    myNet.feedForward(inputVals);

    vector<double> targetVals;
    myNet.backProp(targetVals);

    vector<double> resultVals;
    myNet.getResults(resultVals);
    return 0;
}
