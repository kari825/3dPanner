///
//  BasicSOFA.cpp
//  BasicSOFA
//
//  Created by Allen Lee on 2020-07-09.
//  Copyright © 2020 meoWorkshop. All rights reserved.
//

#include <iostream>
#include <algorithm>
#include "BasicSOFA.hpp"
#include "BasicSOFAPriv.hpp"
#include <JuceHeader.h>

int main() {
    std::cout << "Hello, libBasicSOFA!" << std::endl;
    return 0;
}

inline int roundedKey(double value)
{
    return static_cast<int>(std::round(value * 100.0)); 
}

inline int roundedRadiusKey(double radius)
{
    return static_cast<int>(std::round(radius * 100.0 ));
}

inline int roundedPhiKey(double phi)
{
    return static_cast<int>(std::round(phi * 100.0));
}

inline int roundedThetaKey(double theta)
{
    int key = static_cast<int>(std::round(theta * 100.0));
    if (key > 18000)
        key -= 36000;
    return key;
}

// argMin以上で最小の、argMapのキーの値を返す
inline int findClosestMapKey(int argMin, std::unordered_map<int, size_t> argMap)
{
    const int MAX_KEY = 100000000;
    const int NO_KEY_IN_MAP = -10000000;
    int result = MAX_KEY;
    int keyMax = NO_KEY_IN_MAP;
    for (auto itr = argMap.begin(); itr != argMap.end(); itr++)
    {
        auto key = itr->first;
        // キーが引数の値以上の場合、最終結果が引数の値以上で最小になるように、resultを更新する
        if (key >= argMin)
        {
            if (key < result)
            {
                result = key;
            }
        }
        // キーの最大値を記録しておく
        if (key > keyMax)
        {
            keyMax = key;
        }
    }

    // 引数の値よりも大きい値がない場合、キーの最大値を返す
    if (result == MAX_KEY)
    {
        result = keyMax;
    }

    return result;
}

namespace BasicSOFA
{
    void BasicSOFA::HelloWorld(const char* s)
    {
        BasicSOFAPriv* theObj = new BasicSOFAPriv;
        theObj->HelloWorldPriv(s);
        delete theObj;
    };

    BasicSOFA::BasicSOFA()
    {
        minRadius = 0;
        maxRadius = 0;
        dRadius = 0;
        minPhi = 0;
        maxPhi = 0;
        dPhi = 0;
        minTheta = 0;
        maxTheta = 0;
        dTheta = 0;
        minImpulseDelay = 0;

        fs = 0;
        M = 0;
        N = 0;
        C = 0;
        R = 0;
    }

    bool BasicSOFA::readSOFAFile(std::string filePath)
    {
        resetSOFAData();
        juce::Logger::writeToLog("Trying to open SOFA file:" + juce::String(filePath));
        if (filePath.empty()) {
            juce::Logger::writeToLog("File path is empty, cannot load.");
            return false;
        }

        try
        {
            h5File = H5::H5File(filePath, H5F_ACC_RDONLY);

            H5::Group rootGroup = h5File.openGroup("/");

            // SOFA single dimension parameter size check
            hsize_t dim;

            dim = getSOFASingleDimParameterSize(SOFA_M_STRING);
            if (dim == 0) {
                resetSOFAData();
                return false;
            }
            M = dim;

            dim = getSOFASingleDimParameterSize(SOFA_N_STRING);
            if (dim == 0) {
                resetSOFAData();
                return false;
            }
            N = dim;

            dim = getSOFASingleDimParameterSize(SOFA_R_STRING);
            if (dim == 0) {
                resetSOFAData();
                return false;
            }
            R = dim;

            dim = getSOFASingleDimParameterSize(SOFA_C_STRING);
            if (dim == 0) {
                resetSOFAData();
                return false;
            }
            C = dim;

            std::vector<std::string> possiblePaths = {
                "/Data.SamplingRate",
                "Data.SamplingRate",  
                "/Data/SamplingRate",
                "/SamplingRate"
            };

            bool samplingRateRead = false;

            for (const auto& path : possiblePaths)
            {
                try {
                    H5::DataSet dataSet = h5File.openDataSet(path);
                    H5::DataSpace dataSpace = dataSet.getSpace();

                    int ndims = dataSpace.getSimpleExtentNdims();
                    juce::Logger::writeToLog("Ndims: " + juce::String(ndims));

                    if (ndims != 1)
                    {
                        continue;
                    }

                    hsize_t dims[1] = { 0 };
                    dataSpace.getSimpleExtentDims(dims);
                    if (dims[0] != 1)
                    {
                        continue;
                    }

                    H5::DataType dataType = dataSet.getDataType();
                    H5T_class_t typeClass = dataType.getClass();
                    size_t typeSize = dataType.getSize();

                    if (typeClass == H5T_FLOAT)
                    {
                        if (typeSize == sizeof(double))
                        {
                            double temp = 0.0;
                            H5::DataType memType(H5::PredType::NATIVE_DOUBLE); 
                            dataSet.read(&temp, memType);
                            fs = temp;
                            samplingRateRead = true;

                        }
                        else if (typeSize == sizeof(float))
                        {
                            float temp = 0.0f;
                            H5::DataType memType(H5::PredType::NATIVE_FLOAT);
                            dataSet.read(&temp, memType);
                            fs = static_cast<double>(temp);
                            samplingRateRead = true;

                        }
                    }

                }
                catch (const H5::Exception& e)
                {
                    juce::Logger::writeToLog("Exception reading " + juce::String(path) + ": " + juce::String(e.getDetailMsg().c_str()));
                }
            }

            H5::DataSet dataSet = h5File.openDataSet("Data.SamplingRate");
            H5::DataType dataType = dataSet.getDataType();

            if (dataSet.getStorageSize() > 0) {
                double temp;
                dataSet.read(&temp, dataType);
                fs = temp;
            }
            else {
            }
            



            if (!samplingRateRead)
            {
                try {
                    if (rootGroup.attrExists("SamplingRate")) {
                        H5::Attribute attr = rootGroup.openAttribute("SamplingRate");
                        H5::DataType attrType = attr.getDataType();
                        if (attrType.getClass() == H5T_FLOAT && attrType.getSize() == sizeof(double)) {
                            double temp;
                            attr.read(attrType, &temp);
                            fs = temp;
                            samplingRateRead = true;
                        }
                    }
                    else {
                    }
                }
                catch (const H5::Exception& e) {
                    juce::Logger::writeToLog("Exception reading attribute: " + juce::String(e.getDetailMsg().c_str()));
                }
            }

            if (!samplingRateRead || fs == 0.0) {
                resetSOFAData();
                return false;
            }

            try {
                H5::DataSet irDataset = h5File.openDataSet("Data.IR");
                H5::DataSpace irSpace = irDataset.getSpace();

                hsize_t irDims[3];
                irSpace.getSimpleExtentDims(irDims); // [M][R][N]

                size_t totalSize = irDims[0] * irDims[1] * irDims[2]; // M * R * N
                hrir.resize(totalSize);

                irDataset.read(hrir.data(), H5::PredType::NATIVE_DOUBLE);
            }
            catch (const H5::Exception& e) {
                resetSOFAData();
                return false;
            }

            loadSourcePositions(); 
            bool success = buildCoordinateMap(sourcePositionCoordinates); 
            if (!success)
            {
                resetSOFAData();
                return false;
            }
            radiusMap.clear();
            for (size_t i = 0; i < radiusList.size(); ++i)
            {
                int key = roundedKey(radiusList[i]);
                radiusMap[key] = i;
            }
            dataLoaded = true;
            h5File.close();
            return true;

            dataLoaded = true;
            h5File.close();
            return true;
        }
        catch (H5::Exception& error)
        {
            h5File.close();
            return false;
        }
        catch (H5::DataSetIException& error)
        {
            h5File.close();
            error.printErrorStack();
            return false;
        }
        catch (H5::DataSpaceIException& error)
        {
            h5File.close();
            error.printErrorStack();
            return false;
        }
    }


    
    void BasicSOFA::listAllDatasetNames(const H5::Group& group, const std::string& prefix)
    {

        hsize_t numObjs = group.getNumObjs();
        for (hsize_t i = 0; i < numObjs; ++i)
        {
            std::string name = group.getObjnameByIdx(i);
            int objType = group.getObjTypeByIdx(i);


            if (objType == H5G_GROUP)
            {
                H5::Group subGroup = group.openGroup(name);
                listAllDatasetNames(subGroup, prefix + name + "/");
            }
        }
        
    }

    double getClosest(const std::vector<double>& list, double value)
    {
        auto it = std::min_element(list.begin(), list.end(), [value](double a, double b) {
            return std::abs(a - value) < std::abs(b - value);
            });
        return (it != list.end()) ? *it : value;
    }
    
    const double* BasicSOFA::getHRIR(size_t channel, double theta, double phi, double radius) const noexcept
    {
        if (!dataLoaded || channel >= R) {
            return nullptr;
        }

        if (radiusList.empty() || phiList.empty() || thetaList.empty()) {
            return nullptr;
        }
        
        // 入力値に最も近い値を取ってくる
        double closestRadius = getClosest(radiusList, radius);
        double closestPhi = getClosest(phiList, phi);
        double closestTheta = getClosest(thetaList, theta);

        // radiusから、方位角と仰角にひもづくインデックス番号が入ったMapを取得
        int radiusKey = roundedKey(closestRadius);
        auto radiusIt = radiusMap.find(radiusKey);
        if (radiusIt == radiusMap.end()) {
            return nullptr;
        }
        auto mapIndex = radiusIt->second;
        if (mapIndex >= coordinateMaps.size()) {
            return nullptr;
        }
        const auto& map = coordinateMaps.at(mapIndex); 

        // 仰角のインデックス番号を求める
        int phiKey = roundedPhiKey(closestPhi); 
        auto phiIt = map.phiMap.find(phiKey);
        if (phiIt == map.phiMap.end()) {
            return nullptr;
        }
        size_t phiIndex = phiIt->second;
        if (phiIndex >= map.thetaMaps.size()) {
            return nullptr;
        }

        // 方位角のインデックス番号を求める
        auto& thetaMap = map.thetaMaps.at(phiIndex);
        int thetaKey = roundedThetaKey(closestTheta);
        if (thetaKey > 18000) thetaKey -= 36000;
        int thetaExistingKey = findClosestMapKey(thetaKey, thetaMap);
        auto thetaIt = thetaMap.find(thetaExistingKey);
        if (thetaIt == thetaMap.end()) {
            return nullptr;
        }
        size_t thetaIndex = thetaIt->second;

        // エラーの値がないかチェック
        if (phiIndex >= map.map.size()) {
            return nullptr;
        }
        if (thetaIndex >= map.map[phiIndex].size()) {
            return nullptr;
        }
        if (std::isnan(closestRadius) || std::isnan(closestPhi) || std::isnan(closestTheta))
        {
            return nullptr;
        }

        // 方位角と仰角のインデックス番号から、IRのインデックス番号を求める
        auto irIndex = map.map[phiIndex][thetaIndex];
        if (irIndex * R + channel >= hrir.size()){
            return nullptr;
        }
        return hrir.data() + ((irIndex * R) + channel) * N;
    }

    void BasicSOFA::loadSourcePositions()
    {
        try {

            std::vector<std::string> possiblePaths = {
                "/SourcePosition",
                "SourcePosition",
                "/Data.SourcePosition",
                "Data.SourcePosition",
                "/Data/SourcePosition"
            };

            H5::DataSet dataset;
            bool datasetFound = false;

            for (const auto& path : possiblePaths) {
                try {
                    dataset = h5File.openDataSet(path);
                    datasetFound = true;
                    break;
                }
                catch (const H5::Exception& e) {
                }
            }

            if (!datasetFound) {
                return;
            }

            H5::DataSpace dataspace = dataset.getSpace();
            int ndims = dataspace.getSimpleExtentNdims();
            if (ndims != 2) {
                return;
            }

            hsize_t dims[2];
            dataspace.getSimpleExtentDims(dims);
            hsize_t numPositions = dims[0];
            hsize_t coordCount = dims[1];

            if (coordCount != 3) {
                return;
            }

            std::vector<float> positionData(numPositions * coordCount);
            dataset.read(positionData.data(), H5::PredType::NATIVE_FLOAT);

            auto isClose = [](double a, double b, double eps = 1e-6) {
                return std::abs(a - b) < eps;
            };

            auto existsInList = [&](const std::vector<double>& list, double value) {
                return std::any_of(list.begin(), list.end(), [&](double v) { return isClose(v, value); });
            };

            for (size_t i = 0; i < numPositions; ++i) {
                double azimuth = static_cast<double>(positionData[i * 3 + 0]);
                double elevation = static_cast<double>(positionData[i * 3 + 1]);
                double radius = static_cast<double>(positionData[i * 3 + 2]);

                if (!existsInList(thetaList, azimuth))
                    thetaList.push_back(azimuth);
                if (!existsInList(phiList, elevation))
                    phiList.push_back(elevation);
                if (!existsInList(radiusList, radius))
                    radiusList.push_back(radius);


                if (radiusMap.find(radius) == radiusMap.end()) {
                    coordinateMaps.push_back(SOFACoordinateMap(radius));
                    radiusMap[radius] = coordinateMaps.size() - 1;
                }
            }
            for (size_t i = 0; i < numPositions; ++i) {
                double azimuth = static_cast<double>(positionData[i * 3 + 0]);
                double elevation = static_cast<double>(positionData[i * 3 + 1]);
                double radius = static_cast<double>(positionData[i * 3 + 2]);

                sourcePositionCoordinates.push_back(azimuth);
                sourcePositionCoordinates.push_back(elevation);
                sourcePositionCoordinates.push_back(radius);
            }

        }
        catch (H5::Exception& e) {
            juce::Logger::writeToLog("Exception reading SourcePosition: " + juce::String(e.getDetailMsg().c_str()));
        }
    }


    
    
    void BasicSOFA::resetSOFAData()
    {
        if (hrir.size() != 0)
        {
            hrir.erase(hrir.begin(), hrir.end());
            hrir.shrink_to_fit();
        }
        
        if (radiusList.size() != 0)
        {
            radiusList.erase(radiusList.begin(), radiusList.end());
            radiusList.shrink_to_fit();
        }
        
        if (phiList.size() != 0)
        {
            phiList.erase(phiList.begin(), phiList.end());
            phiList.shrink_to_fit();
        }
        
        if (thetaList.size() != 0)
        {
            thetaList.erase(thetaList.begin(), thetaList.end());
            thetaList.shrink_to_fit();
        }
        
        if (coordinateMaps.size() != 0)
        {
            coordinateMaps.erase(coordinateMaps.begin(), coordinateMaps.end());
            coordinateMaps.shrink_to_fit();
        }
        
        if (radiusMap.size() != 0)
            radiusMap.erase(radiusMap.begin(), radiusMap.end());
        
        minRadius = 0;
        maxRadius = 0;
        dRadius = 0;
        minPhi = 0;
        maxPhi = 0;
        dPhi = 0;
        minTheta = 0;
        maxTheta = 0;
        dTheta = 0;
        
        minImpulseDelay = 0;
        
        fs = 0;
        M = 0;
        N = 0;
        C = 0;
        R = 0;
    }
    
    
    std::vector<double> BasicSOFA::getCoordinatesFromSOFAFile()
    {
        std::vector<double> coordinates;
        
        try
        {
            auto dataSet = h5File.openDataSet(SOFA_SRC_POS_STRING);
            auto dataSpace = dataSet.getSpace();
            
            //  According to the standard, the source and listener positions should have dimensionaltiy [M x C] or [I x C]
            //  The object having M x C dimensionality contains the coordinate data
            hsize_t dims[2];
            dataSpace.getSimpleExtentDims(dims);
            
            if (dims[0] == M)
            {
                coordinates = std::vector<double>(M * C);
                dataSet.read(coordinates.data(), H5::PredType::NATIVE_DOUBLE);
            }
            
            dataSet = h5File.openDataSet(SOFA_LIS_POS_STRING);
            dataSpace = dataSet.getSpace();
            dataSpace.getSimpleExtentDims(dims);
            
            //  If both source and listener objects contain coordinate data, error out and return an empty array
            if (dims[0] == M && coordinates.size() != 0)
            {
                coordinates.erase(coordinates.begin(), coordinates.end());
                coordinates.shrink_to_fit();
                return coordinates;
            }
            
            if (dims[0] == M)
            {
                coordinates = std::vector<double>(M * C);
                dataSet.read(coordinates.data(), H5::PredType::NATIVE_DOUBLE);
            }
            
        }
        catch(...)
        {
            if (coordinates.size() != 0)
            {
                coordinates.erase(coordinates.begin(), coordinates.end());
                coordinates.shrink_to_fit();
            }
            
            return coordinates;
        }
        
        return coordinates;
    }
    
    
    hsize_t BasicSOFA::getSOFASingleDimParameterSize(std::string parameter)
    {
        hsize_t dim;
        
        try
        {
            auto dataSet = h5File.openDataSet(parameter);
            auto dataSpace = dataSet.getSpace();

            dataSpace.getSimpleExtentDims(&dim);
        }
        catch(...)
        {
            return 0;
        }

        return dim;
    }
    
    
    /*
     *  Build tables to map a given radius, theta and phi to its corresponding impulse response
     *
     *  For a given radius, there is a 2D matrix that stores the location of an impulse response for a given theta and phi
     *  Each row of the matrix corresponds to a phi value while each column corresponds to a theta value
     *  For a given phi and theta, there are corresponding std::unordered_map objects that map to the appropriate row and column
     *  This is done to ensure that queries to a particular impulse response is done in constant time - ie: O(1)
     */
    bool BasicSOFA::buildCoordinateMap(const std::vector<double>& coordinates)
    {
        if (coordinates.empty() || (coordinates.size() % C != 0))
            return false;

        auto numBlocks = coordinates.size() / C;

        for (auto block = 0; block < numBlocks; ++block)
        {
            auto thetaIndex = block * C;
            auto phiIndex = block * C +1;
            auto radiusIndex = block * C + 2;

            double theta = coordinates[block * C];
            double phi = coordinates[block * C + 1];

            auto radius = round(coordinates.at(radiusIndex));
           /* juce::Logger::writeToLog("Block " + juce::String(block) +
                ": theta=" + juce::String(coordinates.at(thetaIndex)) +
                ", phi=" + juce::String(coordinates.at(phiIndex)) +
                ", radius=" + juce::String(coordinates.at(radiusIndex)));*/
            if (radius < 0)
                return false;

            int radiusKey = roundedRadiusKey(radius);

            if (radiusMap.find(radiusKey) == radiusMap.end())
            {
                coordinateMaps.push_back(SOFACoordinateMap(radius));
                radiusMap[radiusKey] = coordinateMaps.size() - 1;
                addValueToArray(radius, radiusList);
            }

            auto& map = coordinateMaps[radiusMap[radiusKey]];
           
            int phiKey = roundedPhiKey(coordinates.at(phiIndex));
            if (map.phiMap.find(phiKey) == map.phiMap.end())
            {
                size_t newIndex = map.phiMap.size();
                map.phiMap[phiKey] = newIndex;
                map.thetaMaps.push_back(std::unordered_map<int, size_t>());
                map.map.push_back(std::vector <size_t>());
                addValueToArray(coordinates.at(phiIndex), phiList);
            }

            auto itPhi = map.phiMap.find(phiKey);
            if (itPhi == map.phiMap.end())
            {
                return false;
            }

            auto mapRow = itPhi->second;

            int thetaKey = roundedThetaKey(coordinates.at(thetaIndex));
            if (thetaKey > 18000) thetaKey -= 36000;

            // thetaKeyの存在チェック
            auto& thetaMap = map.thetaMaps[mapRow];
            auto itTheta = thetaMap.find(thetaKey);

            size_t thetaIndexInMap;

            if (itTheta == thetaMap.end())
            {
                thetaIndexInMap = map.map[mapRow].size();  
                thetaMap[thetaKey] = thetaIndexInMap;
                addValueToArray(coordinates.at(thetaIndex), thetaList);
            }
            else
            {
                thetaIndexInMap = itTheta->second;
            }

            if (mapRow >= map.map.size())
            {
                return false;
            }

            auto& mapRowVec = map.map[mapRow];

            if (thetaIndexInMap >= mapRowVec.size())
            {
                mapRowVec.resize(thetaIndexInMap + 1, SIZE_MAX);
            }

            mapRowVec[thetaIndexInMap] = block;
        }

        return true;
    }


    
    
    void BasicSOFA::addValueToArray (const double &x, std::vector<double> &A)
    {
        auto search = std::find(A.begin(), A.end(), x);
        if (search == A.end())
            A.push_back(x);
    }
    
    
    bool BasicSOFA::calculateCoordinateStatisticalData()
    {
        std::sort(radiusList.begin(), radiusList.end());
        std::sort(thetaList.begin(), thetaList.end());
        std::sort(phiList.begin(), phiList.end());
        
        auto delta = round(radiusList[1] - radiusList[0]);
        for (auto i = 2; i < radiusList.size() - 1; ++i)
        {
            if (round(radiusList[i] - radiusList[i - 1]) != delta)
                return false;
        }
        dRadius = delta;
        minRadius = radiusList.at(0);
        maxRadius = radiusList.at(radiusList.size() - 1);
        
        
        delta = abs(round(thetaList[1] - thetaList[0]));
        for (auto i = 2; i < thetaList.size() - 1; ++i)
        {
            if (abs(round(thetaList[i] - thetaList[i - 1])) != delta)
                return false;
        }
        dTheta = delta;
        minTheta = thetaList.at(0);
        maxTheta = thetaList.at(thetaList.size() - 1);
        
        
        delta = abs(round(phiList[1] - phiList[0]));
        for (auto i = 2; i < phiList.size() - 1; ++i)
        {
            if (abs(round(phiList[i] - phiList[i - 1])) != delta)
                return false;
        }
        dPhi = delta;
        minPhi = phiList.at(0);
        maxPhi = phiList.at(phiList.size() - 1);
        
        
        return true;
    }
    
    
    double BasicSOFA::round(const double &x) const
    {
        double temp = 0;
        
        if (x > 0.0)
            temp = x + (epsilon / 2);
        else
            temp = x - (epsilon / 2);
        
        int tempInt = static_cast<int>(temp / epsilon);
        
        return tempInt * epsilon;
    }
    
    
    /*
     *  Go through all of the impulse responses and find where the impulse peak happens and track the index
     *  The minimum index is the 'earliest' impulse delay
     *  Of course, you should not simply use this value when truncating your HRIRs
     *  The starting point of the truncated HRIR should be well less than minImpulseDelay
     *  minImpulseDelay / 2 is a good place to start
     */
    bool BasicSOFA::findMinImpulseDelay()
    {
        if (M == 0)
            return false;
        
        minImpulseDelay = N;
        
        for (auto i = 0; i < M; ++i)
        {
            double hrirMax = 0.0;
            size_t maxLocation = 0;
            
            for (auto j = i * N; j < (i * N) + N; ++j)
            {
                if (hrirMax < abs(hrir[j]))
                {
                    hrirMax = abs(hrir[j]);
                    maxLocation = j;
                }
            }
            
            if (maxLocation < minImpulseDelay)
                minImpulseDelay = maxLocation;
        }
        
        return true;
    }
    
}



void BasicSOFAPriv::HelloWorldPriv (const char * s)
{
    std::cout << s << std::endl;
};


