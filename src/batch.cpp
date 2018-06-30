/*
batch.cpp -- command line interface to Instant Meshes

This file is part of the implementation of

Instant Field-Aligned Meshes
Wenzel Jakob, Daniele Panozzo, Marco Tarini, and Olga Sorkine-Hornung
In ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2015)

All rights reserved. Use of this source code is governed by a
BSD-style license that can be found in the LICENSE.txt file.
*/

#include "batch.h"
#include "meshio.h"
#include "hierarchy.h"
#include "optimizer.h"
#include "dedge.h"
#include "subdivide.h"
#include "timer.h"

void batch_process(char *input, char *output, Float scale) {

	MultiResolutionHierarchy mRes;
	Optimizer *mOptimizer;

	mRes.Two_rosy_flag = true;

	mOptimizer = new Optimizer(mRes);

	Timer<> timer;
	timer.beginStage("data pre-processing");
	mRes.load(input);

	if (scale < 0)
	{
		float tmpScale = mRes.scale();
		float nearest = roundf(tmpScale * 10000) / 10000;
		mRes.setScale(nearest);
	}
	else
		mRes.setScale(scale);

	mRes.build();

	timer.beginStage("rosy optimization");
	
	mOptimizer->setExtrinsic(false);

	mOptimizer->setAlignment(true);
	mOptimizer->setRandomization(true);
	mOptimizer->setHierarchy(true);

	mOptimizer->setOptimizeOrientations(true);
	mOptimizer->notify();
	mOptimizer->wait();

	mRes.detectOrientationSingularitiesTri();

	timer.endStage();
	mRes.sta.timings.push_back(timer.value());

	timer.beginStage("posy optimization");

	mOptimizer->setOptimizePositions(true);
	mOptimizer->notify();
	mOptimizer->wait();

	mRes.detectPositionSingularitiesTri();

	timer.endStage();
	mRes.sta.timings.push_back(timer.value());

	mOptimizer->shutdown();

	timer.beginStage("mesh extraction");

	mRes.re_color = true;
	mRes.splitting = true;
	mRes.decomposes = true;
	mRes.doublets = true;
	mRes.triangles = true;

	if (!mRes.tetMesh()) {
		mRes.re_color = true;
		mRes.splitting = true;
		mRes.decomposes = true;
		mRes.doublets = true;
		mRes.triangles = true;

		mRes.meshExtraction2D();
	}
	else {
		mRes.re_color = true;
		mRes.splitting = true;
		mRes.meshExtraction3D();
	}

	timer.endStage();
	mRes.sta.timings.push_back(timer.value());

	mRes.convert2Poly();
	mRes.labelMesh(false);
	
	mRes.alignMesh();

	mRes.stitchMeshing();

	mRes.exportResult(output);
}
