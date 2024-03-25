#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int extraMemoryAllocated;

void *Alloc(size_t sz)
{
    extraMemoryAllocated += sz;
    size_t* ret = malloc(sizeof(size_t) + sz);
    *ret = sz;
    printf("Extra memory allocated, size: %ld\n", sz);
    return &ret[1];
}

void DeAlloc(void* ptr)
{
    size_t* pSz = (size_t*)ptr - 1;
    extraMemoryAllocated -= *pSz;
    printf("Extra memory deallocated, size: %ld\n", *pSz);
    free((size_t*)ptr - 1);
}

size_t Size(void* ptr)
{
    return ((size_t*)ptr)[-1];
}

void heapify(int arr[], int n, int i) {

    // Largest becomes the root
    int largest = i;
    //Left and then right children are created
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    // If the left child is larger than the root largest
    if (left < n && arr[left] > arr[largest]) {
        largest = left;
    }

    // If the right child is larger than the root largest
    if (right < n && arr[right] > arr[largest]) {
        largest = right;
    }

    // If largest is not root
    if (largest != i) {
        // temp variable swap
        int temp = arr[i];
        arr[i] = arr[largest];
        arr[largest] = temp;

        // Recursively call heapify for the current roots tree
        heapify(arr, n, largest);
    }
}
// implements heap sort
// extraMemoryAllocated counts bytes of memory allocated
void heapSort(int arr[], int start, int end) {

    // Due to how the base program calls heapsort with the beginning and end index
    // I create n here to represent the total number of elements that need sorting in the array
    int n = end + 1;

    // Use heapify to construct the heap and set root to largest
    for (int i = n / 2 - 1; i >= 0; i--) {
        heapify(arr, n, i);
    }

    // Take largest element and place it at root of heap over and over
    for (int i = n - 1; i > 0; i--) {
        // Another temp value swap for the heapSort
        int temp = arr[0];
        arr[0] = arr[i];
        arr[i] = temp;

        // Call heapify to fix the heap on the reduced heap size
        heapify(arr, i, 0);
    }
}

//This will merge the two pData minor arrays created
void merge(int* pData, int l, int m, int r) {
    int i, j, k;
    //This gets the elments in the left and right minor arrays
    int n1 = m - l + 1;
    int n2 = r - m;

    // This allocates memory for two temporary arrays to store the data
    int *L = (int *)Alloc(n1 * sizeof(int)), *R = (int *)Alloc(n2 * sizeof(int));

    // Copy data to the new L and R minor arrays that just got allocated
    for (i = 0; i < n1; i++)
        L[i] = pData[l + i];
    for (j = 0; j < n2; j++)
        R[j] = pData[m + 1+ j];

    // i is the index of first minor array j is the second and k is the index of the newly merged array
    i = 0; j = 0; k = l;
    //Merging of these temporary arrays that are created
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            pData[k] = L[i];
            i++;
        } else {
            pData[k] = R[j];
            j++;
        }
        k++;
    }

    //Copy the elements of L that are left
    while (i < n1) {
        pData[k] = L[i];
        i++;
        k++;
    }

    //Copy the elements of R that are left
    while (j < n2) {
        pData[k] = R[j];
        j++;
        k++;
    }

    //Using our personal DeAlloc function deallocate the memory of the L and R minor arrays created.
    DeAlloc(L);
    DeAlloc(R);
}

// implement merge sort
// extraMemoryAllocated counts bytes of extra memory allocated
void mergeSort(int pData[], int l, int r)
{
    if (l < r) {
        //Here the midpoint is found and this helps divide the array in two
        int m = l+(r-l)/2;

        // Recursive call to sort the first and then second half of the array
        mergeSort(pData, l, m);
        mergeSort(pData, m+1, r);

        //The two halves are sorted by calling them into merge which takes the two sorted halves and merges them together with the two temp arrays in merge
        merge(pData, l, m, r);
    }
}

// implement insertion sort
// extraMemoryAllocated counts bytes of memory allocated
void insertionSort(int* pData, int n)
{
    //Iterate over the array
    int i, key, j;
    for (i = 1; i < n; i++) {
        //Key for the current elment that is to be inserted
        key = pData[i];
        //Start comparing with elements just before the current one which is why i starts at 1 and not 0
        j = i - 1;

        //Elements that are greater than key are moved one forward from their previous position as j moves backwards
        while (j >= 0 && pData[j] > key) {
            pData[j + 1] = pData[j];
            j = j - 1;
        }
        //The key is placed into the proper position
        //Continue doing this with the loop until the array is sorted
        pData[j + 1] = key;
    }
}

// implement bubble sort
// extraMemoryAllocated counts bytes of extra memory allocated
void bubbleSort(int* pData, int n)
{
    //Outer loop that goes through all but the last element
    int i, j;
    // This variable checks to see if any swaps happen
    int swapped;
    for (i = 0; i < n-1; i++) {
        // start swapped at 0 for false
        swapped = 0;

        // This inner loop will compare and potentially swap things
        //As the array area of the array gets smaller keep grabbing the largest element and swapping it
        for (j = 0; j < n-i-1; j++)
            if (pData[j] > pData[j+1]) {
                //Variable temp swap between the j index and what's in front of it
                int temp = pData[j];
                pData[j] = pData[j+1];
                pData[j+1] = temp;
                // swapped variable tells us a swap has occured for the if statement below
                swapped = 1;
            }

        // If no swaps happen the array is sorted already and this is best case
        if (swapped == 0)
            break;
    }
}

// implement selection sort
// extraMemoryAllocated counts bytes of extra memory allocated
void selectionSort(int* pData, int n)
{
    int i, j, minIdx;
    for (i = 0; i < n-1; i++) {
        minIdx = i;
        for (j = i+1; j < n; j++)
            //Loop through the array looking for the current min index
            if (pData[j] < pData[minIdx])
                minIdx = j;

        // temp variable swap of the min element and the current elements position in the array
        int temp = pData[minIdx];
        pData[minIdx] = pData[i];
        pData[i] = temp;
    }

}

// parses input file to an integer array
int parseData(char *inputFileName, int **ppData)
{
    FILE* inFile = fopen(inputFileName, "r");
    int dataSz = 0;
    *ppData = NULL;

    //Was file opened correctly
    if (inFile)
    {
        //Read size of the data
        fscanf(inFile, "%d\n", &dataSz);
        //Using Alloc to allocate memory for the data
        *ppData = (int *)Alloc(sizeof(int) * dataSz);

        // If the memory is allocated incorrectly close out the file and return -1 to indicate failure
        if (!*ppData) {
            //This occurs in case of memory not being allocated correctly with a -1 as the data returned
            fprintf(stderr, "Failed to allocate memory correctly %s\n", inputFileName);
            fclose(inFile);
            return -1;
        }

        //Loop through the file to store each int
        for (int i = 0; i < dataSz; ++i) {
            fscanf(inFile, "%d", (*ppData) + i);
        }

        //Finish storing that data and read throug hthe file then close it out
        fclose(inFile);
    }
    else {
        //This occurs in case of file not opening and willl return this error messagge with a -1 as the data returned
        fprintf(stderr, "Failed to open the file called %s\n", inputFileName);
        return -1;
    }

    return dataSz; // Return the number of data elements read.
}

// prints first and last 100 items in the data array
void printArray(int pData[], int dataSz)
{
	int i, sz = dataSz - 100;
	printf("\tData:\n\t");
	for (i=0;i<100;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\t");

	for (i=sz;i<dataSz;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\n");
}

int main(void)
{
	clock_t start, end;
	int i;
    double cpu_time_used;
	char* fileNames[] = {"input1.txt", "input2.txt", "input3.txt"};

	for (i=0;i<3;++i)
	{
		int *pDataSrc, *pDataCopy;
		int dataSz = parseData(fileNames[i], &pDataSrc);

		if (dataSz <= 0)
			continue;

		pDataCopy = (int *)Alloc(sizeof(int)*dataSz);

		printf("---------------------------\n");
		printf("Dataset Size : %d\n",dataSz);
		printf("---------------------------\n");

		printf("Selection Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		selectionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

		printf("Insertion Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		insertionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

		printf("Bubble Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		bubbleSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

		printf("Merge Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		mergeSort(pDataCopy, 0, dataSz - 1);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

                printf("Heap Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		heapSort(pDataCopy, 0, dataSz - 1);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

		DeAlloc(pDataCopy);
		DeAlloc(pDataSrc);
	}

}
