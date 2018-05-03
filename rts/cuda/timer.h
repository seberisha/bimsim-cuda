static cudaEvent_t tStartEvent;
static cudaEvent_t tStopEvent;

static void gpuStartTimer()
{
	//set up timing events
	cudaEventCreate(&tStartEvent);
	cudaEventCreate(&tStopEvent);
	cudaEventRecord(tStartEvent, 0);
}

static float gpuStopTimer()
{
	cudaEventRecord(tStopEvent, 0);
	cudaEventSynchronize(tStopEvent);
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, tStartEvent, tStopEvent);
	cudaEventDestroy(tStartEvent);
	cudaEventDestroy(tStopEvent);
	return elapsedTime;
}
