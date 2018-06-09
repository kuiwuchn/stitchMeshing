#pragma once
#include "hierarchy.h"
#include <thread>

class Optimizer {
public:
    Optimizer(MultiResolutionHierarchy &mRes);

    void notify() { mCond.notify_all(); }
    void setOptimizeOrientations(bool value);
    void setOptimizePositions(bool value);
    void setAlignment(bool alignment) { mAlignment = alignment; }
    void setRandomization(bool randomization) { mRandomization = randomization; }
    void setExtrinsic(bool extrinsic) { mExtrinsic = extrinsic; }
    void setHierarchy(bool hierarchy) { mHierarchy = hierarchy; }
    void setMaxIterations(uint32_t it) { mMaxIterations = it; }

    bool hierarchy() const { return mHierarchy; }
    bool randomization() const { return mRandomization; }
    bool alignment() const { return mAlignment; }
    bool extrinsic() const { return mExtrinsic; }
    uint32_t maxIterations() const { return mMaxIterations; }

    bool active() const {
        return mRunning && (mOptimizeOrientations || mOptimizePositions);
    }

    void run(); 
	void wait();
    void shutdown();

    void save(Serializer &serializer) const;
    void load(Serializer &serializer);

private:
    MultiResolutionHierarchy &mRes;
    std::thread mThread;
    std::condition_variable_any mCond;
    bool mRunning;
    bool mOptimizeOrientations;
    bool mOptimizePositions;
    bool mAlignment;
    bool mRandomization;
    bool mExtrinsic;
    bool mHierarchy;
    uint32_t mLevel;
    uint32_t mLevelIterations;
    uint32_t mMaxIterations;
};
