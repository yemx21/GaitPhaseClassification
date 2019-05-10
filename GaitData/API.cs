using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace GPC
{
    public static class NativeAPI
    {
        /// <summary>
        /// create a dataset object
        /// </summary>
        /// <param name="channelcount">number of observation channels or number of features per frame</param>
        /// <param name="phasecount">number of phases in a cycle</param>
        /// <returns>address to the dataset object</returns>
        [DllImport("NativeAPI.dll", EntryPoint = "CreateDataSet", CallingConvention = CallingConvention.Cdecl)]
        public static extern IntPtr CreateDataSet(uint channelcount, uint phasecount);

        /// <summary>
        /// destory resource of the dataset object
        /// </summary>
        /// <param name="dataset">address to the dataset object</param>
        [DllImport("NativeAPI.dll", EntryPoint = "DestoryDataSet", CallingConvention = CallingConvention.Cdecl)]
        public static extern void DestoryDataSet(IntPtr dataset);

        /// <summary>
        /// save the dataset object as file
        /// </summary>
        /// <param name="dataset">address to the dataset object</param>
        [DllImport("NativeAPI.dll", EntryPoint = "SaveDataSet", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Unicode)]
        [return: MarshalAs(UnmanagedType.I1)]
        public static extern bool SaveDataSet(IntPtr dataset, string path);

        /// <summary>
        /// create a dataset object from file
        /// </summary>
        /// <param name="channelcount">number of observation channels or number of features per frame</param>
        /// <param name="phasecount">number of phases in a cycle</param>
        /// <returns>address to the dataset object</returns>
        [DllImport("NativeAPI.dll", EntryPoint = "LoadDataSet", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Unicode)]
        public static extern IntPtr LoadDataSet(string path, uint channelcount, uint phasecount);

        /// <summary>
        /// create an empty cycle
        /// </summary>
        /// <param name="complete">true if cycle is a standard full cycle</param>
        /// <returns>address to the cycle object</returns>
        [DllImport("NativeAPI.dll", EntryPoint = "CreateCycle", CallingConvention = CallingConvention.Cdecl)]
        public static extern IntPtr CreateCycle(IntPtr dataset, bool complete);

        /// <summary>
        /// add the cycle object to the dataset object
        /// </summary>
        /// <param name="dataset">address to the dataset object</param>
        /// <param name="cycle">address to the cycle object</param>
        /// <returns></returns>
        [DllImport("NativeAPI.dll", EntryPoint = "AddCycle", CallingConvention = CallingConvention.Cdecl)]
        [return: MarshalAs(UnmanagedType.I1)]
        public static extern bool AddCycle(IntPtr dataset, IntPtr cycle);

        /// <summary>
        /// fill frame into cycle (frames should be standaradized)
        /// </summary>
        /// <param name="cycle">address to the cycle object</param>
        /// <param name="time">timestamp of the frame</param>
        /// <param name="phase">phase of the frame</param>
        /// <param name="vals">observations(features) of the frame</param>
        /// <param name="valcount">number of observations(features)</param>
        /// <returns></returns>
        [DllImport("NativeAPI.dll", EntryPoint = "AddFrame", CallingConvention = CallingConvention.Cdecl)]
        [return: MarshalAs(UnmanagedType.I1)]
        public static extern bool AddFrame(IntPtr cycle, float time, int phase, IntPtr vals, uint valcount);

        /// <summary>
        /// perform gait pattern clustering and averaging
        /// </summary>
        /// <param name="dataset">address to the dataset object</param>
        /// <param name="edgeratio">continuous edge ratio (see Eq.6)</param>
        /// <returns></returns>
        [DllImport("NativeAPI.dll", EntryPoint = "Clustering", CallingConvention = CallingConvention.Cdecl)]
        [return: MarshalAs(UnmanagedType.I1)]
        public static extern bool Clustering(IntPtr dataset, float edgeratio=0.1f);

        /// <summary>
        /// perform optimized feature extraction
        /// </summary>
        /// <param name="dataset">address to the dataset object</param>
        /// <param name="desiredfeaturecount">desired number of resulting features</param>
        /// <param name="duplicate_neigborsize">number of adjacent standardized frames that are used for determining whether feature pairs are duplicated</param>
        /// <returns>address to the featurepairs object</returns>
        [DllImport("NativeAPI.dll", EntryPoint = "ExtractFeaturePairs", CallingConvention = CallingConvention.Cdecl)]
        public static extern IntPtr ExtractFeaturePairs(IntPtr dataset, int desiredfeaturecount, int duplicate_neigborsize);

        /// <summary>
        /// compute and save extracted features
        /// </summary>
        /// <param name="dataset">address to the dataset object</param>
        /// <param name="featurepairs">address to the featurepairs object</param>
        /// <param name="outpath">path to the output file</param>
        /// <param name="standardizecycleframecount">the fixed length that is used to standardize gait pattern (see denote "L")</param>
        /// <returns></returns>
        [DllImport("NativeAPI.dll", EntryPoint = "ComputeFeatures", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Unicode)]
        [return: MarshalAs(UnmanagedType.I1)]
        public static extern bool ComputeFeatures(IntPtr dataset, IntPtr featurepairs, string outpath, uint standardizecycleframecount = 100);
    }
}
