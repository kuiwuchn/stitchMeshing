#include "optimizer.h"

Optimizer::Optimizer(MultiResolutionHierarchy &mRes)
	: mRes(mRes), mRunning(false), mOptimizeOrientations(false),
	mOptimizePositions(false), mAlignment(true), mRandomization(true), mExtrinsic(true),
	mHierarchy(true), mLevel(0), mLevelIterations(0), mMaxIterations(20) {
	mThread = std::thread(&Optimizer::run, this);
}

void Optimizer::setOptimizeOrientations(bool value) {
	std::lock_guard<ordered_lock> lock(mRes.mutex());
	if (mRes.tetMesh()) mMaxIterations = 200; else mMaxIterations = 20;
	mOptimizeOrientations = value;
	mLevel = mRes.levels() - 2;
	mLevelIterations = 0;
}

void Optimizer::setOptimizePositions(bool value) {
	std::lock_guard<ordered_lock> lock(mRes.mutex());
	if (mRes.tetMesh()) mMaxIterations = 200; else mMaxIterations = 20;
	mOptimizePositions = value;
	mLevel = mRes.levels() - 2;
	mLevelIterations = 0;
}

void Optimizer::run() {
	mRunning = true;

	while (true) {
		std::lock_guard<ordered_lock> lock(mRes.mutex());
		while (mRunning && (mRes.levels() == 0 || (!mOptimizePositions && !mOptimizeOrientations)))
			mCond.wait(mRes.mutex());

		if (!mHierarchy)
			mLevel = 0;

		if (mOptimizeOrientations) {
			if (mRes.tetMesh())
				mRes.smoothOrientationsTet(mLevel, mAlignment, mRandomization);
			else
				mRes.smoothOrientationsTri(mLevel, mAlignment, mRandomization, mExtrinsic);
		}

		if (mOptimizePositions) {
			if (mRes.tetMesh())
				mRes.smoothPositionsTet(mLevel, mAlignment, mRandomization);
			else {
				mRes.smoothPositionsTri(mLevel, mAlignment, mRandomization, mExtrinsic);
			}
		}

		mLevelIterations++;

		if (mLevelIterations >= mMaxIterations) {
			mLevelIterations = 0;
			if (mLevel == 0) {
				mOptimizeOrientations = false;
				mOptimizePositions = false;
				notify();
				continue;
			}
			if (mHierarchy) {
				mLevel--;
				if (mOptimizeOrientations)
					mRes.prolongOrientations(mLevel);
				if (mOptimizePositions)
					mRes.prolongPositions(mLevel);
			}
		}

		if (!mRunning)
			break;
	}
}
void Optimizer::wait() {
	std::lock_guard<ordered_lock> lock(mRes.mutex());
	while (mRunning && (mOptimizePositions || mOptimizeOrientations))
		mCond.wait(mRes.mutex());
}
void Optimizer::shutdown() {
	mRunning = false;
	notify();
	mThread.join();
}

void Optimizer::save(Serializer &serializer) const {
	serializer.push("optimizer");
	serializer.set("alignment", mAlignment);
	serializer.set("randomization", mRandomization);
	serializer.set("hierarchy", mHierarchy);
	serializer.set("level", mLevel);
	serializer.set("maxIterations", mMaxIterations);
	serializer.pop();
}

void Optimizer::load(Serializer &serializer) {
	serializer.push("optimizer");
	serializer.get("alignment", mAlignment);
	serializer.get("randomization", mRandomization);
	serializer.get("hierarchy", mHierarchy);
	serializer.get("level", mLevel);
	serializer.get("maxIterations", mMaxIterations);
	serializer.pop();
	mOptimizeOrientations = mOptimizePositions = false;
	mLevelIterations = 0;
}
