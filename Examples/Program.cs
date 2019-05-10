using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using GPC;


namespace Examples
{
    class Program
    {
        static void Main(string[] args)
        {
           
        }

        static IEnumerable<object> GetStandardizedCycles(object datafile)
        {
            //unimplemented
            return new List<object>();
        }

        static IEnumerable<object> GetNormalCycles(object datafile)
        {
            //unimplemented
            return new List<object>();
        }

        static bool IsComplete(object cycle)
        {
            //unimplemented
            return true;
        }

        static IEnumerable<ValueTuple<float, int, float[]>> GetFrames(object cycle)
        {
            //unimplemented
            return new List<ValueTuple<float, int, float[]>>();
        }

        /// <summary>
        /// load standardized dataset
        /// </summary>
        /// <param name="datafiles">some files</param>
        /// <param name="Z">number of features per frame</param>
        /// <param name="K">number of gait phases</param>
        static IntPtr LoadStandardizedDataset(List<object> datafiles, uint Z = 12, uint K = 9)
        {
            // create an empty dataset object
            var dataset = NativeAPI.CreateDataSet(Z, K);

            foreach (var any_datafile in datafiles)
            {
                foreach (var any_cycle in GetStandardizedCycles(any_datafile))
                {
                    // create an empty cycle object
                    var cycle = NativeAPI.CreateCycle(dataset, IsComplete(any_cycle));

                    foreach (var any_frame in GetFrames(cycle))
                    {
                        //fill frame to the cycle object
                        NativeAPI.AddFrame(cycle, any_frame.Item1, any_frame.Item2, Marshal.UnsafeAddrOfPinnedArrayElement(any_frame.Item3, 0), (uint)any_frame.Item3.Length);
                    }

                    //add the cycle object to dataset
                    NativeAPI.AddCycle(dataset, cycle);
                }
            }
            return dataset;
        }


        /// <summary>
        /// load training dataset
        /// </summary>
        /// <param name="datafiles">some files</param>
        /// <param name="Z">number of features per frame</param>
        /// <param name="K">number of gait phases</param>
        static IntPtr LoadNormalDataset(List<object> datafiles, uint Z=12, uint K=9)
        {
            // create an empty dataset object
            var dataset = NativeAPI.CreateDataSet(Z, K);

            foreach (var any_datafile in datafiles)
            {
                foreach (var any_cycle in GetNormalCycles(any_datafile))
                {
                    // create an empty cycle object
                    var cycle = NativeAPI.CreateCycle(dataset, IsComplete(any_cycle));

                    foreach (var any_frame in GetFrames(cycle))
                    {
                        //fill frame to the cycle object
                        NativeAPI.AddFrame(cycle, any_frame.Item1, any_frame.Item2, Marshal.UnsafeAddrOfPinnedArrayElement(any_frame.Item3, 0), (uint)any_frame.Item3.Length);
                    }

                    //add the cycle object to dataset
                    NativeAPI.AddCycle(dataset, cycle);
                }
            }
            return dataset;
        }


        static void Example1(List<object> trainingdatafiles, List<object> testingdatafiles, string outtrainpath, string outtestpath)
        {
            //number of frames in standardized gait cycle curves
            uint L = 100;
            //continuous edge ratio
            float tau = 0.1f;
            //number of gait phases
            uint K = 9;
            //number of features per frame
            uint Z = 12;
            //total desired number of feature pairs
            int Ohm = 420;
            //number of adjacent standardized frames that are used for determining whether feature pairs are duplicated
            int t = 6;


            //Step1: fill in standardized cycles from training set
            var dataset = LoadStandardizedDataset(trainingdatafiles, Z, K);

            //Step2: perform gait pattern clustering and averaging
            NativeAPI.Clustering(dataset, tau);

            //Step3: perform feature pair extraction
            var featurepairs = NativeAPI.ExtractFeaturePairs(dataset, Ohm, t);
            NativeAPI.DestoryDataSet(dataset);//release resources if not used anymore


            //Step5: compute and export features of the training dataset to file
            // choose one from the following two lines according to data scale
            var training_dataset = LoadNormalDataset(trainingdatafiles, Z, K);
            //var training_dataset = LoadStandardizedDataset(trainingdatafiles, Z, K);
            NativeAPI.ComputeFeatures(training_dataset, featurepairs, outtrainpath, L);
            NativeAPI.DestoryDataSet(training_dataset);//release resources if not used anymore

            //Step6: compute and export features of the testing dataset to file
            var testing_dataset = LoadNormalDataset(testingdatafiles, Z, K);
            NativeAPI.ComputeFeatures(testing_dataset, featurepairs, outtestpath, L);
            NativeAPI.DestoryDataSet(testing_dataset);//release resources if not used anymore
        }

        static void Example2(List<object> trainingdatafiles, List<object> testingdatafiles, string outtrainpath, string outtestpath)
        {
            //number of gait phases
            uint K = 9;
            //number of features per frame
            uint Z = 12;

            //fill data in a dataset
            var dataset = LoadNormalDataset(trainingdatafiles, Z, K);

            //save the dataset to somewhere
            NativeAPI.SaveDataSet(dataset, @"drive:\someplace");
            NativeAPI.DestoryDataSet(dataset);//release resources if not used anymore


            //reload the dataset from somewhere
            var dataset1 = NativeAPI.LoadDataSet(@"drive:\someplace", Z, K);
            NativeAPI.DestoryDataSet(dataset1);//release resources if not used anymore
        }
    }
}
