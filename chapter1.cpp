#include<vector>

template <typename T>
std::vector<T> operator*(const std::vector<T>& vec, T factor) {
    std::vector<T> result = vec;
    for (auto& value : result) {
        value *= factor;
    }
    return result;
}
// Overload * operator to scale the values (scale * vector)
template <typename T>
std::vector<T> operator*(T factor, const std::vector<T>& vec) {
    std::vector<T> result = vec;
    for (auto& value : result) {
        value *= factor;
    }
    return result;
}
template <typename T>
std::vector<T> operator+(const std::vector<T>& vec1, const std::vector<T>& vec2) {
    if (vec1.size() != vec2.size()) {
        throw std::runtime_error("Vector sizes must be the same for addition.");
    }

    std::vector<T> result;
    result.reserve(vec1.size());

    for (size_t i = 0; i < vec1.size(); ++i) {
        result.push_back(vec1[i] + vec2[i]);
    }

    return result;
}


class Sys{};

float getTime(Sys sys){
    return 1.0f;
}

std::vector<float> deriveEval(Sys sys, std::vector<float> x0, float t){
    return x0;
}
std::vector<float> getState(Sys sys){
    return std::vector<float>{0.0f};
}
void setState(Sys sys, std::vector<float> xh,  float t){

}
void eurlerStep(Sys sys, float h){
    float t = getTime(sys);
    std::vector<float> x0, deltaX;

    t = getTime(sys);
    x0 = getState(sys);
    deltaX =  deriveEval(sys, x0, t); 
    setState(sys,x0 + h * deltaX, t + h );
}   