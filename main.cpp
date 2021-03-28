#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <chrono>

#define TYPE int
#define LENGTH 1000000

using namespace std;

template<typename T>
void Merge(T arr[], long long start, long long end)
{
    long long z, x, y, mid;
    vector<T> temp(end - start + 1);
    mid = (start + end) / 2;
    z = 0;
    x = start;
    y = mid + 1;

    while (x <= mid && y <= end)
    {
        if (arr[x] < arr[y])
        {
            temp[z] = arr[x];
            ++x, ++z;
        }
        else
        {
            temp[z] = arr[y];
            ++y, ++z;
        }
    }

    while (x <= mid)
    {
        temp[z] = arr[x];
        ++x, ++z;
    }

    while (y <= end)
    {
        temp[z] = arr[y];
        ++y, ++z;
    }

    for (long long i = start; i <= end; ++i)
    {
        arr[i] = temp[i - start];
    }

}

template<typename T>
void MergeSort(T arr[], long long start, long long end)
{
    if (start < end)
    {
        long long mid = (start + end) / 2;
        MergeSort(arr, start, mid);
        MergeSort(arr, mid + 1, end);
        Merge(arr, start, end);
    }

}

template<typename T>
long long Partition(T arr[], long long start, long long end) {
	long long pivot = arr[end];
	long long temp;
	long long i = start;

	for (long long j = start; j < end; ++j)
	{
		if (arr[j] <= pivot)
		{
			temp = arr[j];
			arr[j] = arr[i];
			arr[i] = temp;
			i++;
		}
	}

	arr[end] = arr[i];
	arr[i] = pivot;

	return i;
}

template<typename T>
void QuickSort(T arr[], long long start, long long end)
{
	if (start < end)
	{
		long long q = Partition(arr, start, end);
		QuickSort(arr, start, q - 1);
		QuickSort(arr, q + 1, end);
	}
}

template<typename T>
void InsertionSort(T arr[], long long end)
{
    long long i, key, j;
    for (i = 1; i < end; i++) {
        key = arr[i];
        j = i - 1;

        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

template<typename T>
void MaxHeapify(T arr[], long long heapSize, long long index) {
	long long left = (index + 1) * 2 - 1;
	long long right = (index + 1) * 2;
	long long largest = 0;

	if (left < heapSize && arr[left] > arr[index])
		largest = left;
	else
		largest = index;

	if (right < heapSize && arr[right] > arr[largest])
		largest = right;

	if (largest != index)
	{
		long long temp = arr[index];
		arr[index] = arr[largest];
		arr[largest] = temp;

		MaxHeapify(arr, heapSize, largest);
	}
}

template<typename T>
void HeapSort(T arr[], long long end) {
	long long heapSize = end;

	for (long long p = (heapSize - 1) / 2; p >= 0; --p)
		MaxHeapify(arr, heapSize, p);

	for (long long i = end - 1; i > 0; --i)
	{
		long long temp = arr[i];
		arr[i] = arr[0];
		arr[0] = temp;

		--heapSize;
		MaxHeapify(arr, heapSize, 0);
	}
}

template<typename T>
void IntroSort(T arr[], long long end)
{
	long long partitionSize = Partition(arr, 0, end);

	if (partitionSize < 16)
	{
		InsertionSort(arr, end);
	}
	else if (partitionSize >(2 * log(end)))
	{
		HeapSort(arr, end);
	}
	else
	{
		QuickSort(arr, 0, end - 1);
	}
}

template<typename T>
void PrintArray(T arr[], long long n)
{
    for (long long i = 0; i < n; ++i)
        cout << arr[i] << " ";
    cout << endl << endl;
}

template<typename T>
void RandomizeArray(T arr[], long long n)
{
    T MIN=(-1000000);
    T MAX=1000000;
    srand(time(NULL));
    for(long long i=0;i<n;i++)
    {
        arr[i]=MIN + (double)rand() / RAND_MAX * (MAX - MIN);
    }
}

template<typename T>
void FillArrayInRisingOrder(T arr[], long long p, long long n)  //p - up to which point array will be filled in order
{
    if(p>n)
    {
        cerr << "Error with FillArrayInOrder" << endl;
        exit(-1);
    }
    for(int i = 0; i<p; i++)
    {
        arr[i]=i;
    }
}

template<typename T>
void FillArrayInLoweringOrder(T arr[], long long p, long long n)  //p - up to which point array will be filled in order
{
    if(p>n)
    {
        cerr << "Error with FillArrayInOrder" << endl;
        exit(-1);
    }
    for(int i = 0; i<p; i++)
    {
        arr[i]=-i;
    }
}

template<typename T>
void ReverseArray(T arr[], long long end)  //Works wrong with odd numbers but it is not needed in used lengths of array, plus it wasn't used
{
    for(int i = 0; i <= end/2; i++)
    {
        swap(arr[i], arr[end-i]);
    }
}

template<typename T>
bool CheckArray (T arr[], long long n)
{
    T prev=arr[0];
    for(long long i=1; i<n ;i++)
    {
        T curr=arr[i];
        if(prev>curr)
        {
            return 0;
        }
        prev=curr;
    }
    return 1;
}

template<typename T>
void MakeArray(T arr[], long long end)
{
    //RandomizeArray(arr, end);
    FillArrayInLoweringOrder(arr, end, end);
    //FillArrayInRisingOrder(arr, 0.25*end, end);
}

int main()
{
    using micro_s = std::chrono::microseconds;
    //using milli_s = std::chrono::milliseconds;
    //using seconds = std::chrono::seconds;

    long long n = LENGTH;
    TYPE *arr = new TYPE[LENGTH];


    //cout << "Array Before Sorting: " << endl;
    //PrintArray(arr, n);

    std::chrono::time_point<std::chrono::high_resolution_clock> tsum[100][3][2];

    for(int i=0; i<100; i++)
    {
        MakeArray(arr, n);

        auto t1 = std::chrono::high_resolution_clock::now();

        MergeSort(arr, 0, n-1);

        auto t2 = std::chrono::high_resolution_clock::now();

        tsum[i][0][0]=t1;
        tsum[i][0][1]=t2;

        MakeArray(arr, n);

        t1 = std::chrono::high_resolution_clock::now();

        QuickSort(arr, 0, n - 1);

        t2 = std::chrono::high_resolution_clock::now();

        tsum[i][1][0]=t1;
        tsum[i][1][1]=t2;

        MakeArray(arr, n);

        t1 = std::chrono::high_resolution_clock::now();

        IntroSort(arr, n);

        t2 = std::chrono::high_resolution_clock::now();

        tsum[i][2][0]=t1;
        tsum[i][2][1]=t2;

    }

    //cout << "Array After Sorting: " << endl;
    //PrintArray(arr, n);

    //if(!CheckArray(arr, n))
    //{
    //    cerr << "Sorting failed" << endl;
    //}
    //else
    //{
    //    cout << "Sorting succesfull!" << endl;
    //}

    long long TM = 0, TQ = 0, TI = 0;

    for(int i = 0; i < 100; i++)
    {
        TM=TM+std::chrono::duration_cast<micro_s>( tsum[i][0][1] - tsum[i][0][0] ).count();
        TQ=TQ+std::chrono::duration_cast<micro_s>( tsum[i][1][1] - tsum[i][1][0] ).count();
        TI=TI+std::chrono::duration_cast<micro_s>( tsum[i][2][1] - tsum[i][2][0] ).count();
    }
    cout << "Scalanie: " << (double)TM/1000000 << " s" << endl;
    cout << "QuickSort: " << (double)TQ/1000000 << " s" << endl;
    cout << "IntroSort: " << (double)TI/1000000 << " s" << endl;

    delete arr;
    arr=NULL;

    return 0;
}
