#include "pch.h"
#include <aef/operators/IOperator.h>
#include <unordered_map>

namespace {
    std::unordered_map<std::string, aef::operators::IOperator::pfnOperatorMaker> &opTypeMap = registerDefaults();

    std::unordered_map<std::string, aef::operators::IOperator::pfnOperatorMaker> &registerDefaults() {
        auto* map = new std::unordered_map<std::string, aef::operators::IOperator::pfnOperatorMaker>;
        

        return *map;
    }
}

void aef::operators::IOperator::registerOperatorType(std::string name, pfnOperatorMaker ctor) {
    opTypeMap[name] = ctor;
}

aef::operators::IOperator* aef::operators::IOperator::makeOperatorOfType(std::string name) {
    if (!opTypeMap.contains(name)) {
        return nullptr;
    }
    auto pfn = opTypeMap[name];
    return pfn();
}
