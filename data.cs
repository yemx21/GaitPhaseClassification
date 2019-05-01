using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Media3D;
using System.Xml.Serialization;

namespace GPC
{
    public enum MarkerPosition
    {
        Invalid,

        Left_Shoulder,
        Right_Shoulder,

        Left_Anterior_Superior_Iliac_Spine,
        Right_Anterior_Superior_Iliac_Spine,

        Left_Posterior_Superior_Iliac_Spine,
        Right_Posterior_Superior_Iliac_Spine,

        Left_Hip,
        Right_Hip,

        Left_Fema,
        Right_Fema,

        Left_Knee,
        Right_Knee,

        Left_Tibial,
        Right_Tibial,

        Left_Ankle,
        Right_Ankle,

        Left_Toe,
        Right_Toe,

        Left_Heel,
        Right_Heel,
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct Point3d
    {
        public float X;
        public float Y;
        public float Z;
        public Point3d(float x, float y, float z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        public Point3d(Point3d pt)
        {
            X = pt.X;
            Y = pt.Y;
            Z = pt.Z;
        }

        public Point3d(double x, double y, double z)
        {
            X = (float)x;
            Y = (float)y;
            Z = (float)z;
        }

        public static Point3d Empty;
        public static Point3d NaN;

        static Point3d()
        {
            Empty = new Point3d(0, 0, 0);
            NaN = new Point3d(double.NaN, double.NaN, double.NaN);
        }

        public Vector3D ToVector3D()
        {
            return new Vector3D(X, Y, Z);
        }

        public bool HasNaN()
        {
            return double.IsNaN(X) || double.IsNaN(Y) || double.IsNaN(Z);
        }
    }


    [StructLayout(LayoutKind.Sequential)]
    public struct Point2d
    {
        public float X;
        public float Y;
        public Point2d(float x, float y)
        {
            X = x;
            Y = y;
        }
    }

    [StructLayout(LayoutKind.Sequential)]
    [Serializable]
    public class MarkerData
    {
        public float X;
        public float Y;
        public float Z;
        public MarkerPosition Pos;

        /* OX, OY are discard */
        public float OX;
        public float OY;

        public MarkerData() { X = float.NaN; Y = float.NaN; Z = float.NaN; Pos = MarkerPosition.Invalid; OX = float.NaN; OY = float.NaN; }
        public MarkerData(float x, float y, float z, float ox, float oy, MarkerPosition pos) { X = x; Y = y; Z = z; OX = ox; OY = oy; Pos = pos; }
        public MarkerData(Point3d xyz, Point2d oxy, MarkerPosition pos) { X = xyz.X; Y = xyz.Y; Z = xyz.Z; OX = oxy.X; OY = oxy.Y; Pos = pos; }
        public MarkerData(Point3d xyz, MarkerPosition pos) { X = xyz.X; Y = xyz.Y; Z = xyz.Z; OX = float.NaN; OY = float.NaN; Pos = pos; }


        public static bool IsNaNOrNull(MarkerData md)
        {
            if (md == null) return true;
            return double.IsNaN(md.X) || double.IsNaN(md.Y) || double.IsNaN(md.Z);

        }

        public Point3d GetPosition()
        {
            return new Point3d(X, Y, Z);
        }

        public void Mute()
        {
            X = float.NaN;
            Y = float.NaN;
            Z = float.NaN;
        }

    }

    [Serializable]
    public class GaitFrame
    {
        public List<MarkerData> Markers = new List<MarkerData>();

        public float Time { get; set; } = -1.0f;

        public GaitFrame() { }
        public GaitFrame(GaitFrame other)
        {
            Time = other.Time;
            Markers = other.Markers.ConvertAll(mk => new MarkerData(mk.X, mk.Y, mk.Z, mk.OX, mk.OY, mk.Pos));
        }

        public GaitFrame(FrameData fd)
        {
            Time = fd.TimeIndex;
            Markers = fd.Markers.ConvertAll(mk => new MarkerData(mk.X, mk.Y, mk.Z, mk.OX, mk.OY, mk.Pos));
        }
    }

    [Serializable]
    public class GaitFloor
    {
        /*A,B,C means normal vector (x,y,z) of the floor plane, D is the y coordinate refers to camera */
        public double A { get; set; }
        public double B { get; set; }
        public double C { get; set; }
        public double D { get; set; }

        public GaitFloor()
        {

        }

        public GaitFloor(double a, double b, double c, double d)
        {
            A = a;
            B = b;
            C = c;
            D = d;
        }

        public double DistanceFrom(Point3d pos3, bool isSigned = false)
		{
			return isSigned? (A* pos3.X + B* pos3.Y + C* pos3.Z + D) * (B< 0.0 ? 1 : -1) : Math.Abs(A* pos3.X + B* pos3.Y + C* pos3.Z + D);
        }   
    }

    public enum GaitEvent : int
    {
        Invalid = 0,
        InitialContact = 1,
        FootFlat,
        Toe1Off,
        Midstance,
        HeelOff,
        Balance,
        ToeOff,
        Midswing,
        Foot1Only,
    }

    public enum GaitPhase : int
    {
        Invalid = 0,
        LoadingResponse1 = 1,
        LoadingResponse2,
        Midstance,
        Terminalstance1,
        Terminalstance2,
        Preswing,
        Initialswing,
        Midswing,
        Terminalswing
    }

    [Serializable]
    public class GaitEventShot
    {
        public GaitEvent Event { get; set; } = GaitEvent.Invalid;

        public float Time { get; set; } = -1.0f;

        public GaitEventShot()
        {

        }
    }

    public enum Sagittal
    {
        Left,
        Right
    }

    [Serializable]
    public class GaitData
    {
        public List<GaitFrame> Frames = new List<GaitFrame>();

        public List<GaitEventShot> Events = new List<GaitEventShot>();

        public Sagittal Sagittal { get; set; } = Sagittal.Left;

        public GaitFloor Floor { get; set; } = null;

        public int SubjectId { get; set; } = -1;

        public Guid ExperimentId { get; set; } = Guid.Empty;

    }
}

