#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import copy
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import itertools
import shapely
import shapely.ops
import time
from shapely.geometry import Polygon
import json

#==========================================================================================================================================================

class xGeoUtils:
  @staticmethod
  def rotate2D(SrcXY, Angle):
    c, s = np.cos(Angle), np.sin(Angle)
    RotMat = np.matrix(((c, -s), (s, c)), dtype=np.float64)
    DstXY  = RotMat * np.asmatrix(SrcXY).T
    Resulr = np.asarray(DstXY).reshape(2)
    return Resulr

  @staticmethod
  def normalizeVector(SrcXY):
    Norm = np.linalg.norm(SrcXY)
    return SrcXY/Norm

  @staticmethod
  def lineEquationCoeffsFromTwoPoints(A, B): 
    #ax + by + c = 0
    a = B[1] - A[1]
    b = A[0] - B[0]
    c = A[1]*B[0] - A[0]*B[1]
    return (a, b, c)

  @staticmethod
  def calcPointToLineDistance(PointXY, LineCoeffs):
    a, b, c = LineCoeffs
    x, y    = PointXY
    Distance = abs(a*x + b*y + c) / math.sqrt(a**2 + b**2)
    return Distance

  @staticmethod
  def findIntersection(LineCoeffsA, LineCoeffsB):
    try:
      a = np.array([[LineCoeffsA[0], LineCoeffsA[1]], [LineCoeffsB[0], LineCoeffsB[1]]])
      b = np.array([-LineCoeffsA[2], -LineCoeffsB[2]])
      R = np.linalg.solve(a, b)
      return R
    except np.linalg.LinAlgError:
      return []

  @staticmethod
  def crossProduct2D(PointXY, VectorBegXY, VectorEndXY):
    CrossProduct = (VectorEndXY[0] - VectorBegXY[0])*(PointXY[1] - VectorBegXY[1]) - (VectorEndXY[1] - VectorBegXY[1])*(PointXY[0] - VectorBegXY[0])
    return CrossProduct

  @staticmethod
  def isWithinTriangle(PointXY, Triangle):
    VertexA, VertexB, VertexC = Triangle
    CrossProductAB = xGeoUtils.crossProduct2D(PointXY, VertexA, VertexB)
    CrossProductBC = xGeoUtils.crossProduct2D(PointXY, VertexB, VertexC)
    CrossProductCA = xGeoUtils.crossProduct2D(PointXY, VertexC, VertexA)
    if(CrossProductAB>=0 and CrossProductBC>=0 and CrossProductCA>=0): return True
    if(CrossProductAB< 0 and CrossProductBC< 0 and CrossProductCA< 0): return True
    return False

  @staticmethod
  def calcAngleBetwenVectors(LeftVector, RightVector):
    Angle = math.atan2(np.linalg.det([LeftVector, RightVector]), np.dot(LeftVector, RightVector))
    return Angle

#==========================================================================================================================================================

class xCamera:
  def __init__(self, ID):
    #params
    self.ID          = ID
    self.FocalLength = None #f
    self.SensorWidth = None #d
    self.Resolution  = None
    self.PixelPitch  = None
    self.AngleHalf   = None
    self.AngleFull   = None
    #ccords
    self.CoordXY     = None
    self.Rotation    = None
    self.Distance    = None
    self.CameraVersorXY    = None
    #rays
    self.Rays          = None
    self.Intersections = {}
    self.Regions       = []
    self.RegDistances  = []
    self.RegErrors     = []
    return

  def __str__(self):
    return "CAM{id} f={f} CoordXY={x},{y} Rot={r}".format(id=self.ID, f=self.FocalLength, x=self.CoordXY[0], y=self.CoordXY[1], r=math.degrees(self.Rotation))

  def __repr__(self):
    return "<<{}>>".format(str(self))

  def setParams(self, FocalLength, SensorWidth, Resolution):
    self.FocalLength = FocalLength
    self.SensorWidth = SensorWidth
    self.Resolution  = Resolution
    self.PixelPitch  = self.SensorWidth / self.Resolution
    self.AngleHalf   = math.atan(self.SensorWidth / (2*self.FocalLength))
    self.AngleFull   = 2*self.AngleHalf        
    return

  def setCoords(self, CoordXY, Rotation):
    self.CoordXY      = np.asarray(CoordXY)
    self.Rotation     = Rotation         

  def calculateDerrived(self):
    self.CameraVersorXY   = xGeoUtils.rotate2D((1,0), self.Rotation)
    self.SensorVersorXY   = xGeoUtils.rotate2D(self.CameraVersorXY, math.radians(90))
    self.RealSensorCntrXY = self.CoordXY - self.CameraVersorXY * self.FocalLength
    self.RealSensorBegXY  = self.RealSensorCntrXY - self.SensorVersorXY * self.SensorWidth / 2
    self.RealSensorEndXY  = self.RealSensorCntrXY + self.SensorVersorXY * self.SensorWidth / 2
    self.FakeSensorCntrXY = self.CoordXY + self.CameraVersorXY * self.FocalLength
    self.FakeSensorBegXY  = self.FakeSensorCntrXY - self.SensorVersorXY * self.SensorWidth / 2
    self.FakeSensorEndXY  = self.FakeSensorCntrXY + self.SensorVersorXY * self.SensorWidth / 2

    self.OpticalAxisCoeffs = xGeoUtils.lineEquationCoeffsFromTwoPoints(self.CoordXY, self.CoordXY + self.CameraVersorXY)
    self.CameraPlaneCoeffs = xGeoUtils.lineEquationCoeffsFromTwoPoints(self.CoordXY, self.CoordXY + self.SensorVersorXY)

    self.DistanceToCamera  = np.linalg.norm(self.CoordXY) #math.sqrt(self.CoordXY[0]**2 + self.CoordXY[1]**2)     
    self.DistanceToPlane   = self.calcDistanceToCameraPlane([0, 0])

    self.RealPixelEdgeCoordsXY = [self.RealSensorBegXY + i * self.SensorVersorXY * self.PixelPitch for i in range(0, self.Resolution + 1)]
    self.FakePixelEdgeCoordsXY = [self.FakeSensorBegXY + i * self.SensorVersorXY * self.PixelPitch for i in range(0, self.Resolution + 1)]

    self.RaysVersor = [xGeoUtils.normalizeVector(self.FakePixelEdgeCoordsXY[i] - self.CoordXY) for i in range(0, self.Resolution + 1)]
    self.RaysCoeffs = [xGeoUtils.lineEquationCoeffsFromTwoPoints(self.CoordXY, self.FakePixelEdgeCoordsXY[i]) for i in range(0, self.Resolution + 1)]
    return

  def calcDistanceToOpticalCenter(self, PointXY):
    return np.linalg.norm(self.CoordXY - PointXY)

  def calcDistanceToCameraPlane(self, PointXY):
    return xGeoUtils.calcPointToLineDistance(PointXY, self.CameraPlaneCoeffs)

  def calcDistanceToCameraOpticalAxis(self, PointXY):
    return xGeoUtils.calcPointToLineDistance(PointXY, self.OpticalAxisCoeffs)

  def isBehindCamera(self, PointXY):
    CrossProduct = xGeoUtils.crossProduct2D(PointXY, self.FakeSensorBegXY, self.FakeSensorEndXY)
    return CrossProduct > 0

  def isWithinFieldOfView(self, PointXY):
    FieldOfViewTriangle = self.getFieldOfViewTriangle()
    Result = xGeoUtils.isWithinTriangle(PointXY, FieldOfViewTriangle)
    return Result

  def calculateIntersections(self, OtherCamera):
      Intersections = {}
      for i, MyRay in enumerate(self.RaysCoeffs):
          MyRayIntersections = {}
          for j, OtherRay in enumerate(OtherCamera.RaysCoeffs):                
              Intersection = xGeoUtils.findIntersection(MyRay, OtherRay)
              if(len(Intersection)):
                  if(not self.isBehindCamera(Intersection)):
                      MyRayIntersections[j] = Intersection
          if MyRayIntersections:
              Intersections[i] = MyRayIntersections
      self.Intersections[OtherCamera.ID] = Intersections
      return

  def calculateRegions(self):
    self.Regions = []
    for CamId in self.Intersections.keys():
      for MyRayId in self.Intersections[CamId].keys():
        if(not MyRayId+1 in self.Intersections[CamId]): continue
        for OtherRayId in self.Intersections[CamId][MyRayId].keys():
          if(not OtherRayId+1 in self.Intersections[CamId][MyRayId  ]): continue
          if(not OtherRayId+1 in self.Intersections[CamId][MyRayId+1]): continue     
          if(not OtherRayId   in self.Intersections[CamId][MyRayId+1]): continue     
          
          RegionCoords = ((self.Intersections[CamId][MyRayId    ][OtherRayId    ]),
                          (self.Intersections[CamId][MyRayId    ][OtherRayId + 1]),
                          (self.Intersections[CamId][MyRayId + 1][OtherRayId + 1]),
                          (self.Intersections[CamId][MyRayId + 1][OtherRayId    ]))
          self.Regions.append(RegionCoords)
    return

  def calculateDistancesAndErrors(self):
    self.Distances = []
    self.Errors    = []
    for i, Region in enumerate(self.Regions):
      Distances = [xGeoUtils.calcPointToLineDistance(VertexXY, self.CameraPlaneCoeffs) for VertexXY in Region]
      MinDistance = np.min(Distances)
      MaxDistance = np.max(Distances)
      self.RegDistances.append(np.average    (Distances))
      self.RegErrors   .append(MaxDistance - MinDistance)

  def getFieldOfViewTriangle(self):
    DefaultLength = self.DistanceToPlane*20
    TargetL = self.CoordXY + (np.array(xGeoUtils.rotate2D(self.CameraVersorXY,  self.AngleHalf)) * DefaultLength)
    TargetR = self.CoordXY + (np.array(xGeoUtils.rotate2D(self.CameraVersorXY, -self.AngleHalf)) * DefaultLength)
    return (self.CoordXY, TargetL, TargetR)

  def calcPixelFieldOfViewVersorsFlex(self, PointXY):
    if(self.isBehindCamera(PointXY)): return None
    RayLineCoeffs    = xGeoUtils.lineEquationCoeffsFromTwoPoints(PointXY, self.CoordXY)
    SensorLineCoeffs = xGeoUtils.lineEquationCoeffsFromTwoPoints(self.FakeSensorBegXY, self.FakeSensorEndXY)
    IntersectionXY   = xGeoUtils.findIntersection(RayLineCoeffs, SensorLineCoeffs)
    if(not len(IntersectionXY)): return None
    DistanceFakeSensorBeg = np.linalg.norm(self.FakeSensorBegXY - IntersectionXY)
    DistanceFakeSensorEnd = np.linalg.norm(self.FakeSensorEndXY - IntersectionXY)
    if(DistanceFakeSensorBeg>self.SensorWidth or DistanceFakeSensorEnd>self.SensorWidth): return None
    LeftXY  = IntersectionXY - self.SensorVersorXY * (self.PixelPitch / 2)
    RightXY = IntersectionXY + self.SensorVersorXY * (self.PixelPitch / 2)
    LeftVersor  = xGeoUtils.normalizeVector(LeftXY  - self.CoordXY)
    RightVersor = xGeoUtils.normalizeVector(RightXY - self.CoordXY)
    return (LeftVersor, RightVersor)

  def calcPixelFieldOfViewTriangleFlex(self, PointXY):
    Versors = self.calcPixelFieldOfViewVersorsFlex(PointXY)
    if(Versors == None): return None
    LeftVersor, RightVersor = Versors
    DefaultLength = self.DistanceToPlane*20
    Triangle = (self.CoordXY, self.CoordXY + LeftVersor*DefaultLength, self.CoordXY + RightVersor*DefaultLength)
    return Triangle

  def calcPixelFieldOfViewAngle(self, PointXY):
    Versors = self.calcPixelFieldOfViewVersorsFlex(PointXY)
    if(Versors == None): return None
    LeftVersor, RightVersor = Versors
    Angle = math.atan2(np.linalg.det([LeftVersor,RightVersor]),np.dot(LeftVersor,RightVersor))
    return Angle

  def getExactPixelTriangle(self, PointXY):
    if(self.isBehindCamera(PointXY)): return None
    RayLineCoeffs    = xGeoUtils.lineEquationCoeffsFromTwoPoints(PointXY, self.CoordXY)
    SensorLineCoeffs = xGeoUtils.lineEquationCoeffsFromTwoPoints(self.FakeSensorBegXY, self.FakeSensorEndXY)
    Intersection     = xGeoUtils.findIntersection(RayLineCoeffs, SensorLineCoeffs)
    if(not len(Intersection)): return None
    DistanceFakeSensorBeg = np.linalg.norm(self.FakeSensorBegXY - Intersection)
    DistanceFakeSensorEnd = np.linalg.norm(self.FakeSensorEndXY - Intersection)
    if(DistanceFakeSensorBeg>self.SensorWidth or DistanceFakeSensorEnd>self.SensorWidth): return None
    IdxA = math.floor(DistanceFakeSensorBeg / self.PixelPitch)
    IdxB = math.ceil (DistanceFakeSensorBeg / self.PixelPitch)
    IdxA = np.clip(IdxA, 0, self.Resolution)
    IdxB = np.clip(IdxB, 0, self.Resolution)
    if(IdxA == IdxB):
      if  (IdxA == 0              ): IdxB = 1
      elif(IdxB == self.Resolution): IdxA = self.Resolution-1
      else: IdxB = IdxB + 1

    DefaultLength = self.DistanceToPlane*20
    Triangle = (self.CoordXY, self.CoordXY+self.RaysVersor[IdxA]*DefaultLength, self.CoordXY+self.RaysVersor[IdxB]*DefaultLength)
    return Triangle

  def drawCamera(self, Axes, Color):   
    DefaultLength = self.DistanceToPlane*2
    #optical center
    plt.scatter(self.CoordXY[0], self.CoordXY[1], s = DefaultLength/200, color = Color)
    #central line
    TargetC = self.CoordXY + self.CameraVersorXY * DefaultLength
    plt.plot([self.CoordXY[0], TargetC[0]], [self.CoordXY[1], TargetC[1]], color = Color, linewidth=1.5, linestyle="--")
    #AOV lines
    AngleOfViewLength = DefaultLength / math.cos(self.AngleHalf)
    TargetL = self.CoordXY + (np.array(xGeoUtils.rotate2D(self.CameraVersorXY, self.AngleHalf)) * AngleOfViewLength)
    plt.plot([self.CoordXY[0], TargetL[0]], [self.CoordXY[1], TargetL[1]], color = Color, linewidth=1.5, linestyle="-")
    TargetR = self.CoordXY + (np.array(xGeoUtils.rotate2D(self.CameraVersorXY, -self.AngleHalf)) * AngleOfViewLength)
    plt.plot([self.CoordXY[0], TargetR[0]], [self.CoordXY[1], TargetR[1]], color = Color, linewidth=1.5, linestyle="-")
    #camera body
    Unit  = -self.FocalLength*2
    Poly0 = self.CoordXY
    Poly1 = self.CoordXY + (np.array(xGeoUtils.rotate2D(self.CameraVersorXY, -self.AngleHalf)) * Unit)
    Poly2 = self.CoordXY + (np.array(xGeoUtils.rotate2D(self.CameraVersorXY,  self.AngleHalf)) * Unit)
    PolyBody = plt.Polygon([Poly0, Poly1, Poly2], edgecolor = Color, facecolor="none")
    Axes.add_patch(PolyBody)
    #camera sensor
    plt.plot([self.RealSensorBegXY[0], self.RealSensorEndXY[0]], [self.RealSensorBegXY[1], self.RealSensorEndXY[1]], color = "black")
    return

  def drawRays(self, Axes, Color):
      DefaultLength = self.DistanceToPlane*2
      plt.plot([self.FakeSensorBegXY[0], self.FakeSensorEndXY[0]], [self.FakeSensorBegXY[1], self.FakeSensorEndXY[1]], color = "black")
      for i in range(0, self.Resolution + 1):
          FakePixelEdgeCoordXY = self.FakePixelEdgeCoordsXY[i]
          RayVector = FakePixelEdgeCoordXY - self.CoordXY
          RayVector = RayVector / np.linalg.norm(RayVector)
          RayVector = RayVector * DefaultLength
          RayPoint  = FakePixelEdgeCoordXY + RayVector
          plt.plot([FakePixelEdgeCoordXY[0], RayPoint[0]], [FakePixelEdgeCoordXY[1], RayPoint[1]], linewidth=0.5, color = Color)
      return

  def drawRayIntersections(self, Axes, Color):
    for IntersectionsCam in self.Intersections.values():
      for IntersectionsMyRay in IntersectionsCam.values():
        for IntersectionOtherRay in IntersectionsMyRay.values():        
          plt.scatter(IntersectionOtherRay[0], IntersectionOtherRay[1], s = self.DistanceToPlane/200, color = Color)
    return

  def drawRegions(self, Axes):
    cmap = plt.get_cmap("tab20")
    for i, Region in enumerate(self.Regions):
      PolyBody = plt.Polygon(Region, edgecolor="none", facecolor=cmap(i%cmap.N))
      Axes.add_patch(PolyBody)
    return

  def drawRegionErrors(self, Axes):
    if not self.RegDistances: return
    cmap = plt.get_cmap("jet")
    ColorMin = 0
    ColorMax = cmap.N
    Selector = self.RegDistances < self.DistanceToPlane*2
    SelErr   = np.asarray(self.RegErrors)[Selector]
    ErrorMin = np.min(SelErr)
    ErrorMax = np.max(SelErr)
    RegionColorIndex = ColorMin + ((self.RegErrors-ErrorMin)/(ErrorMax-ErrorMin))*ColorMax
    for i, Region in enumerate(self.Regions):
      if(self.RegDistances[i] < self.DistanceToPlane*2):
        PolyBody = plt.Polygon(Region, color=cmap(int(RegionColorIndex[i])))
        Axes.add_patch(PolyBody)
    return

#==========================================================================================================================================================

class xArrangement:
  def __init__(self):
    self.FocalLength = None #f
    self.SensorWidth = None #d
    self.Resolution  = None    
    self.ErrorMap    = None    

  def setCameraParams(self, FocalLength, SensorWidth, Resolution):
    self.FocalLength = FocalLength
    self.SensorWidth = SensorWidth
    self.Resolution  = Resolution

  def generateLinear(self, Distance, Baseline, NumCams):
    Cameras  = []
    CoordX   = -Distance 
    CoordY   = Baseline * (NumCams-1) / 2        
    Rotation = math.radians(0)
    for i in range(0, NumCams):
      Camera = xCamera(i)
      Camera.setParams(self.FocalLength, self.SensorWidth, self.Resolution)
      Camera.setCoords((CoordX, CoordY - i*Baseline), Rotation)
      Camera.calculateDerrived()
      Cameras.append(Camera)
    return Cameras
  

  def generateLinearKK(self, Distance, Baseline, NumCams):
    Cameras  = []
    CoordX   = -Distance                     
    CoordY   = Baseline * (NumCams-1) / 2  
    Rotation = math.radians(0)               
    for i in range(0, NumCams):
      Camera = xCamera(i)
      Camera.setParams(self.FocalLength, self.SensorWidth, self.Resolution)
      Camera.setCoords((CoordX+ (np.random.random()-0.5)*Baseline/10, CoordY - i*Baseline + (np.random.random()-0.5)*Baseline/10), Rotation + (np.random.random()-0.5)*0.1)
      Camera.calculateDerrived()
      Cameras.append(Camera)
    return Cameras

  def generateCircular(self, Radius, Angle, NumCams):
    Cameras      = []
    BaseCoordXY  = (-Distance, 0)
    BaseRotation = math.radians(0)
    BegRotation  = -Angle * (NumCams-1) / 2
    for i in range(0, NumCams):
      Rotation = BegRotation + i*Angle
      CoordXY  = xGeoUtils.rotate2D(BaseCoordXY, Rotation)
      Camera = xCamera(i)
      Camera.setParams(self.FocalLength, self.SensorWidth, self.Resolution)
      Camera.setCoords(CoordXY, BaseRotation + Rotation)
      Camera.calculateDerrived()
      Cameras.append(Camera)
    return Cameras

  def generateSportHall(self, Radius, Angle, NumCams):
    Cameras      = []
    #TODO
    return Cameras

#==========================================================================================================================================================

class xEstimator:
  def __init__(self, Cameras):
    self.Cameras   = Cameras
    return

  #polygon method
  def _calcIntersection(self, IsWithinFieldOfView, PointXY, BoxXlim, BoxYlim):  
    Triangles   = [Camera.calcPixelFieldOfViewTriangleFlex(PointXY) for c, Camera in enumerate(self.Cameras) if IsWithinFieldOfView[c]]
    TrianglesOK = [i for i in Triangles if i]
    if(len(TrianglesOK) < 2): return None
    Polygons     = [shapely.geometry.Polygon(Triangle) for Triangle in TrianglesOK]
    Intersection = shapely.geometry.box(BoxXlim[0], BoxYlim[0], BoxXlim[1], BoxYlim[1])
    for Polygon in Polygons: Intersection = Intersection.intersection(Polygon)
    return tuple(Intersection.boundary.coords)

  def calcIntersection(self, PointXY, BoxXlim, BoxYlim):
    IsWithinFieldOfView = [Camera.isWithinFieldOfView(PointXY) for Camera in self.Cameras]
    return self._calcIntersection(IsWithinFieldOfView, PointXY, BoxXlim, BoxYlim)

  def calcDeltaDistances(self, PointXY, BoxXlim, BoxYlim):
    NumOfCameras = len(self.Cameras)
    IsWithinFieldOfView = [Camera.isWithinFieldOfView(PointXY) for Camera in self.Cameras]
    if(len(IsWithinFieldOfView)        <= 1): return np.full((NumOfCameras), np.finfo(np.float64).max, dtype=np.float64)
    if(IsWithinFieldOfView.count(True) <= 1): return np.full((NumOfCameras), np.finfo(np.float64).max, dtype=np.float64)
    IntersectionVertices = self._calcIntersection(IsWithinFieldOfView, PointXY, BoxXlim, BoxYlim)
    DeltaDistances = np.full((NumOfCameras), np.finfo(np.float64).max, dtype=np.float64)
    for c, Camera in enumerate(self.Cameras):
      if(IsWithinFieldOfView[c] and IntersectionVertices ):
        Distances = [ Camera.calcDistanceToCameraPlane(VertexXY) for VertexXY in IntersectionVertices]
        MinDistance = np.min(Distances)
        MaxDistance = np.max(Distances)
        DeltaDistance = MaxDistance - MinDistance
        DeltaDistances[c] = DeltaDistance
    return DeltaDistances

  def __calculateDeltaDistanceMapsPolygonRow(self):
    return

  def calculateDeltaDistanceMapsPolygon(self, LimX, LimY, Precision):
    print("BEG --> calculateDeltaDistanceMaps")
    TimeBeg = time.monotonic()

    NumOfCameras = len(self.Cameras)
    X = np.linspace(LimX[0], LimX[1], Precision)
    Y = np.linspace(LimY[0], LimY[1], Precision)

    RangeX = np.abs(LimX[1] - LimX[0])
    RangeY = np.abs(LimY[1] - LimY[0])
    BoxLimX = (LimX[0] - RangeX/10, LimX[1] + RangeX/10)
    BoxLimY = (LimY[0] - RangeY/10, LimY[1] + RangeY/10)

    DeltaDistanceMaps = [np.full((Precision, Precision), np.finfo(np.float64).max, dtype=np.float64) for i in range(NumOfCameras)]

    for j, y in enumerate(Y): 
      print("{} ".format(j), end="")
      for i, x in enumerate(X):
        PointXY = np.array([x, y])
        DeltaDistances = self.calcDeltaDistances(PointXY, BoxLimX, BoxLimY)
        for c, Camera in enumerate(self.Cameras):
          DeltaDistanceMaps[c][j][i] = DeltaDistances[c]
    print("")

    TimeEnd = time.monotonic()
    print("END --> calculateDeltaDistanceMaps TIME={}".format(TimeEnd - TimeBeg))
    return DeltaDistanceMaps
  
  #SimplifiedForLinear method
  def calculateDeltaDistanceMapsSimplifiedForLinear(self, Xlim, Ylim, Precision):
    print("BEG --> calculateDeltaDistanceMapsSimplifiedForLinear")
    TimeBeg = time.monotonic()

    X = np.linspace(Xlim[0], Xlim[1], Precision)
    Y = np.linspace(Ylim[0], Ylim[1], Precision)

    PixelPitch  = self.Cameras[0].PixelPitch
    FocalLength = self.Cameras[0].FocalLength

    CamPairs = list(itertools.combinations(self.Cameras, 2))
    DeltaZs  = [np.full((Precision, Precision), np.finfo(np.float64).max, dtype=np.float64) for CamPair in range(len(CamPairs))]    

    for p, CamPair in enumerate(CamPairs):
      CameraA, CameraB = CamPair
      Baseline = np.linalg.norm(CameraA.CoordXY - CameraB.CoordXY)
      for j, y in enumerate(Y): 
        for i, x in enumerate(X):
          PointXY = np.array([x, y])          
          if(not CameraA.isWithinFieldOfView(PointXY)): continue
          if(not CameraB.isWithinFieldOfView(PointXY)): continue          
          Z = CameraA.calcDistanceToCameraPlane(PointXY)  
          dZ = 2*abs((PixelPitch * Z**2) / (FocalLength * Baseline - Z * PixelPitch))
          #dZ = abs((PixelPitch * Z**2) / (FocalLength * Baseline - Z * PixelPitch)) 
          DeltaZs[p][j][i] = dZ

    DeltaDistanceMaps = [np.full((Precision, Precision), np.finfo(np.float64).max, dtype=np.float64) for i in range(len(self.Cameras))]
    DeltaCamMaps =      [np.full((Precision, Precision), 255, dtype=np.uint8) for i in range(len(self.Cameras))]

    for c, Camera in enumerate(self.Cameras):
      for p, CamPair in enumerate(CamPairs):
        if(Camera in CamPair):
          print(CamPair[1])
          #DeltaDistanceMaps[c] = np.minimum(DeltaDistanceMaps[c], DeltaZs[p])
          for j, y in enumerate(Y): 
            for i, x in enumerate(X):
              if (DeltaZs[p][i][j] < DeltaDistanceMaps[c][i][j]):
                DeltaDistanceMaps[c][i][j] = DeltaZs[p][i][j]
                if(CamPair[0] == Camera):
                  DeltaCamMaps[c][i][j] = CamPair[1].ID
                else:
                  DeltaCamMaps[c][i][j] = CamPair[0].ID

    TimeEnd = time.monotonic()
    print("END --> calculateDeltaDistanceMapsSimplifiedForLinear TIME={}".format(TimeEnd - TimeBeg))
    return (DeltaDistanceMaps, DeltaCamMaps)

  #SimplifiedForLinear method
  def recalculateDeltaDistanceMapsSimplifiedForLinearAndAllCameras(self, Xlim, Ylim, Precision, DeltaCamMaps):
    print("BEG --> recalculateDeltaDistanceMapsSimplifiedForLinearAndAllCameras")
    TimeBeg = time.monotonic()

    X = np.linspace(Xlim[0], Xlim[1], Precision)
    Y = np.linspace(Ylim[0], Ylim[1], Precision)

    DeltaCamMapsAll =  np.full((Precision, Precision), 0, dtype=np.uint8) 

    for c, Camera in enumerate(self.Cameras):
      print(Camera.ID)
      for j, y in enumerate(Y):
        for i, x in enumerate(X):
          if DeltaCamMaps[c][i][j] == 255 :
            continue
          dist = abs(int(Camera.ID) - int(DeltaCamMaps[c][i][j]) )
          if( dist > DeltaCamMapsAll[i][j] ):
            DeltaCamMapsAll[i][j] = dist   

    DeltaCamMapsAll[DeltaCamMapsAll == 0] = 255

    TimeEnd = time.monotonic()
    print("END --> recalculateDeltaDistanceMapsSimplifiedForLinearAndAllCameras TIME={}".format(TimeEnd - TimeBeg))
    return (DeltaCamMapsAll)

  def calculateDeltaDistanceMapsSimplifiedForCircular(self, Xlim, Ylim, Precision):
    print("BEG --> calculateDeltaDistanceMapsSimplifiedForCircular")
    TimeBeg = time.monotonic()

    X = np.linspace(Xlim[0], Xlim[1], Precision)
    Y = np.linspace(Ylim[0], Ylim[1], Precision)

    PixelPitch  = self.Cameras[0].PixelPitch
    FocalLength = self.Cameras[0].FocalLength

    CamPairs = list(itertools.combinations(self.Cameras, 2))
    #[CurrCam][2ndCam]
    DeltaZs  = [[np.full((Precision, Precision), np.finfo(np.float64).max, dtype=np.float64) for c1 in range(len(self.Cameras))] for c2 in range(len(self.Cameras)) ]  

    for p, CamPair in enumerate(CamPairs):
      CameraA, CameraB = CamPair
      DistanceCamToCam = np.linalg.norm(CameraA.CoordXY - CameraB.CoordXY)
      for j, y in enumerate(Y): 
        for i, x in enumerate(X):
          PointXY = np.array([x, y])          
          if(not CameraA.isWithinFieldOfView(PointXY)): continue
          if(not CameraB.isWithinFieldOfView(PointXY)): continue

          DistancePointToCamA = CameraA.calcDistanceToOpticalCenter(PointXY)
          DistancePointToCamB = CameraB.calcDistanceToOpticalCenter(PointXY)
          DistanceToCamAxisA  = CameraA.calcDistanceToCameraOpticalAxis(PointXY)
          DistanceToCamAxisB  = CameraB.calcDistanceToCameraOpticalAxis(PointXY)

          # cam axi correction

          kat_widzenia_punkt_camA = math.degrees(xGeoUtils.calcAngleBetwenVectors(CameraA.CameraVersorXY, PointXY - CameraA.CoordXY))
          kat_widzenia_punkt_camB = math.degrees(xGeoUtils.calcAngleBetwenVectors(CameraB.CameraVersorXY, PointXY - CameraB.CoordXY))

          kat_widzenia_camB_camA = math.degrees(xGeoUtils.calcAngleBetwenVectors(CameraA.CameraVersorXY, CameraB.CameraVersorXY))
          kat_widzenia_camA_camB = math.degrees(xGeoUtils.calcAngleBetwenVectors(CameraB.CameraVersorXY, CameraA.CameraVersorXY))
          
          if CameraA.calcPixelFieldOfViewAngle(PointXY) == None : continue
          fov_step_camA = math.degrees(CameraA.calcPixelFieldOfViewAngle(PointXY))

          if CameraB.calcPixelFieldOfViewAngle(PointXY) == None : continue
          fov_step_camB = math.degrees(CameraB.calcPixelFieldOfViewAngle(PointXY))

          dif_beta = math.degrees(math.acos((DistancePointToCamA*DistancePointToCamA + DistancePointToCamB*DistancePointToCamB - DistanceCamToCam*DistanceCamToCam)/(2*DistancePointToCamA*DistancePointToCamB))) #% 90

          h_camA = DistancePointToCamA * math.tan(math.radians(fov_step_camA))
          h_camB = DistancePointToCamB * math.tan(math.radians(fov_step_camB))
          

          DeltaZ_CamB_usingCamA = h_camA/math.sin(math.radians(dif_beta)) # omni
          DeltaZ_CamA_usingCamB = h_camB/math.sin(math.radians(dif_beta)) # omni


          DeltaZ_CamB_usingCamA = DeltaZ_CamB_usingCamA * math.cos(math.radians(kat_widzenia_punkt_camB)) # pers
          DeltaZ_CamA_usingCamB = DeltaZ_CamA_usingCamB * math.cos(math.radians(kat_widzenia_punkt_camA)) # persp

  
          IndexA = self.Cameras.index(CameraA)
          IndexB = self.Cameras.index(CameraB)

          DeltaZs[IndexA][IndexB][j][i] = DeltaZ_CamB_usingCamA*2 
          DeltaZs[IndexB][IndexA][j][i] = DeltaZ_CamA_usingCamB*2

    DeltaDistanceMaps = [np.full((Precision, Precision), np.finfo(np.float64).max, dtype=np.float64) for i in range(len(self.Cameras))]
    DeltaCamMaps =  [np.full((Precision, Precision), 255, dtype=np.uint8) for i in range(len(self.Cameras))]
    
    for c1, Camera1 in enumerate(self.Cameras):
      for c2, Camera2 in enumerate(self.Cameras):
        if(c1 != c2):
          for j, y in enumerate(Y):
            for i, x in enumerate(X):
              if  DeltaZs[c1][c2][j][i] < DeltaDistanceMaps[c1][j][i]:
                DeltaDistanceMaps[c1][j][i] =  DeltaZs[c1][c2][j][i]
                DeltaCamMaps[c1][j][i] =  Camera2.ID

    TimeEnd = time.monotonic()
    print("END --> calculateDeltaDistanceMapsSimplifiedForCircular TIME={}".format(TimeEnd - TimeBeg))
    return (DeltaDistanceMaps, DeltaCamMaps)

#==========================================================================================================================================================

class xVisualiser:
  DefaultAxisLabel = "Distance from center of the scene [m]"

  def __init__(self, Estimator):
    self.Estimator = Estimator
    self.Cameras   = Estimator.Cameras
    return

  def setLabels(self, Axes):
    ticks = mpl.ticker.FuncFormatter(lambda x, pos: "{}".format(x/1000))
    Axes.xaxis.set_major_formatter(ticks)
    Axes.yaxis.set_major_formatter(ticks)
    Axes.set_xlabel(self.DefaultAxisLabel)
    Axes.set_ylabel(self.DefaultAxisLabel)
    return

  def drawSystem(self, Axes):    
    Axes.scatter(0, 0, s = 40, color = "black")
    cmap = plt.get_cmap("tab10") if len(self.Cameras)<=10 else plt.get_cmap("tab20")
    for i, Camera in enumerate(self.Cameras):
      CameraColor = cmap(i)
      Camera.drawCamera(Axes, CameraColor)
    self.setLabels(Axes)
    return

  def drawCamDistances(self, Axes, Resolution):    
    Axes.scatter(0, 0, s = 40, color = "black")
    cmap = plt.get_cmap("tab10_r") if len(self.Cameras)<=10 else plt.get_cmap("tab20_r")

    CamPairs = list(itertools.combinations(self.Cameras, 2))

    offset =  Resolution * 0.6


    for p, CamPair in enumerate(CamPairs):
      CameraA, CameraB = CamPair
      if(CameraA.ID != 0): continue
      LinieColor = cmap(abs(CameraA.ID - CameraB.ID))
      plt.annotate('', xy=(CameraA.CoordXY[0]- offset,CameraA.CoordXY[1]), xytext=(CameraB.CoordXY[0]- offset,CameraB.CoordXY[1]), color = "blue", arrowprops=dict(arrowstyle='<->',color = LinieColor,lw=2))
      offset += Resolution * 0.2
    return

  def drawIntersection(self, Axes, IntersectionPoint, Xlim, Ylim):    
    cmap = plt.get_cmap("tab10") if len(self.Cameras)<=10 else plt.get_cmap("tab20")
    for i, Camera in enumerate(self.Cameras):
      CameraColor = cmap(i)      
      Triangle = Camera.calcPixelFieldOfViewTriangleFlex(IntersectionPoint)
      if(Triangle):
        FaceColor = list(CameraColor[:3]) + [0.2]
        PolyBody = plt.Polygon(Triangle, edgecolor=CameraColor, facecolor=FaceColor, linestyle=":")
        Axes.add_patch(PolyBody)
    Intersection = self.Estimator.calcIntersection(IntersectionPoint, Xlim, Ylim)
    if(Intersection):
      PolyBody = plt.Polygon(Intersection, color="black")
      Axes.add_patch(PolyBody)
    return      

  def drawErrorMap(self, Axes, ErrorMap, ErrorMapXlim, ErrorMapYlim):
    cmap = plt.get_cmap("jet")
    ColorMin = 0
    ColorMax = cmap.N
    Selector = ErrorMap != np.finfo(np.float64).max
    SelErr   = np.asarray(ErrorMap)[Selector]
    ErrorMin = np.min(SelErr)
    ErrorMax = np.max(SelErr)
    ColorIndex = ColorMin + ((ErrorMap-ErrorMin)/(ErrorMax-ErrorMin))*ColorMax
    ErrorImg   = np.zeros((ErrorMap.shape[0], ErrorMap.shape[1], 4))
    for y in range(ErrorMap.shape[1]):
      for x in range(ErrorMap.shape[0]):
        if(ColorIndex[y][x] >= ColorMin and ColorIndex[y][x] < ColorMax):
            Color = tuple(cmap(int(ColorIndex[y][x])))
            ErrorImg[y,x,:] = Color
  
    Axes.imshow(ErrorImg, origin="lower", extent=(ErrorMapXlim[0], ErrorMapXlim[1], ErrorMapYlim[0], ErrorMapYlim[1]))#, aspect="auto")  
    norm = mpl.colors.Normalize(vmin=ErrorMin, vmax=ErrorMax)    
    mappable = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = Axes.figure.colorbar(mappable, ax=ax)
    cbar.set_label("Distance error [mm]")
    cbar.ax.locator_params(nbins=10)
    return mappable

  def drawDiffMap(self, Axes, ErrorMap, Mask, ErrorMapXlim, ErrorMapYlim):
    cmap = plt.get_cmap("jet")
    ColorMin = 0
    ColorMax = cmap.N
    Selector = Mask != np.finfo(np.float64).max
    SelErr   = np.asarray(ErrorMap)[Selector]
    ErrorMin = np.min(SelErr)
    ErrorMax = np.max(SelErr)
    ColorIndex = ColorMin + ((ErrorMap-ErrorMin)/(ErrorMax-ErrorMin))*ColorMax
    ErrorImg   = np.zeros((ErrorMap.shape[0], ErrorMap.shape[1], 4))
    for y in range(ErrorMap.shape[1]):
      for x in range(ErrorMap.shape[0]):
        if(ColorIndex[y][x] >= ColorMin and ColorIndex[y][x] <ColorMax):
            Color = tuple(cmap(int(ColorIndex[y][x])))
            ErrorImg[y,x,:] = Color
  
    Axes.imshow(ErrorImg, origin="lower", extent=(ErrorMapXlim[0], ErrorMapXlim[1], ErrorMapYlim[0], ErrorMapYlim[1]))#, aspect="auto")  
    norm = mpl.colors.Normalize(vmin=ErrorMin, vmax=ErrorMax)    
    mappable = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = Axes.figure.colorbar(mappable, ax=ax)
    cbar.set_label("Distance error [mm]")
    cbar.ax.locator_params(nbins=10)
    return

  def drawCamMap(self, Axes, CamMap, ErrorMapXlim, ErrorMapYlim, ColorReverse = False):

    cmap = plt.get_cmap("tab10") if len(self.Cameras)<=10 else plt.get_cmap("tab20_r")
    if ColorReverse :
      cmap = plt.get_cmap("tab10_r") if len(self.Cameras)<=10 else plt.get_cmap("tab20_r")
    ColorMin = 0
    ColorMax = cmap.N
    Selector = CamMap != 255
    SelErr   = np.asarray(CamMap)[Selector]
    ErrorMin = np.min(SelErr)
    ErrorMax = np.max(SelErr)

    print("ErrorMax " +  str(ErrorMax))
    
    ErrorImg   = np.zeros((CamMap.shape[0], CamMap.shape[1], 4))
    
    for i, Camera in enumerate(self.Cameras):
      CameraColor = cmap(i)
      Camera.drawCamera(Axes, CameraColor)
      
      for y in range(CamMap.shape[1]):
        for x in range(CamMap.shape[0]):
          if(CamMap[y][x] != 255):
            Color = tuple(cmap(int(CamMap[y][x])))
            ErrorImg[y,x,:] = Color 
  
    Axes.imshow(ErrorImg, origin="lower", extent=(ErrorMapXlim[0], ErrorMapXlim[1], ErrorMapYlim[0], ErrorMapYlim[1]))#, aspect="auto")  
    norm = mpl.colors.Normalize(vmin=ErrorMin, vmax=ErrorMax)    

    return

#==========================================================================================================================================================


mpl.rcParams.update({'font.size': 20})

f = open("config.json")
jsonparam = json.load(f)
f.close

#camera parameters
FocalLength = jsonparam['FocalLength_mm'] # 35mm
SensorWidth = jsonparam['SensorWidth_mm'];  #36mm
Resolution  = jsonparam['Resolution_px'] #16 pixels total

#scene parameters
Distance  = jsonparam['DisplayDistance_mm'] #1500 mm

json_cir = jsonparam['Circular']
json_lin = jsonparam['Linear']

GenCir = json_cir['Generate']
GenLin = json_lin['Generate']

Angle     =  math.radians(json_cir['CameraAngle_deg']) #math.radians(30)
Baseline  = json_lin['CameraBasline_mm'] #180 Distance * math.sin(Angle)   mm
NumCams   = jsonparam['NumberOfCameras'] #2
MainCamId = jsonparam['MainCamera'] #int(NumCams/2)

#plot parameters
Xlim=(-1.2*Distance, 1.2*Distance)
Ylim=(-1.2*Distance, 1.2*Distance)

#create cameras
Arrangement = xArrangement()
Arrangement.setCameraParams(FocalLength, SensorWidth, Resolution)
CamerasL    = Arrangement.generateLinear(Distance, Baseline, NumCams)
CamerasC    = Arrangement.generateCircular(Distance, Angle, NumCams)
EstimatorL  = xEstimator (CamerasL)
EstimatorC  = xEstimator (CamerasC)
VisualiserL = xVisualiser(EstimatorL)
VisualiserC = xVisualiser(EstimatorC)

DefaultFigSize            = (12,12)
DisplayFigs               = jsonparam['DisplayFigures']
DrawSystemOverview        = jsonparam['DrawSystemOverview']
DrawErrorMap              = jsonparam['DrawErrorMap']
DrawCamMap                = jsonparam['DrawCameraMap']

DrawSimplifiedErrLinear   = jsonparam['DrawSimplified'] & GenLin
DrawSimplifiedErrCircular = jsonparam['DrawSimplified'] & GenCir

ErrorMapResolution = jsonparam['ErrorMapResolution']

if DrawSystemOverview:
  print("DrawSystemOverview...")

  if GenLin:
    fig = plt.figure(figsize=DefaultFigSize)
    ax  = fig.add_subplot(1,1,1)
    ax.set_title("Considered system")
    VisualiserL.drawSystem(ax)
    ax.set(xlim=Xlim, ylim=Ylim)
    ax.grid()
    plt.tight_layout()
    if DisplayFigs: plt.show()
  
  if GenCir:
    fig = plt.figure(figsize=DefaultFigSize)
    ax  = fig.add_subplot(1,1,1)
    ax.set_title("Considered system")
    VisualiserC.drawSystem(ax)
    ax.set(xlim=Xlim, ylim=Ylim)
    ax.grid()
    plt.tight_layout()
    if DisplayFigs: plt.show()    

if DrawErrorMap:
  print("DrawErrorMap...")
  ErrorMapXlim = (-Distance, Distance)
  ErrorMapYlim = (- 1.2 * Distance, 1.2*Distance)
  print("calculateErrorMap...")
  
  if GenLin:
    ErrorMapsPolyL = EstimatorL.calculateDeltaDistanceMapsPolygon(ErrorMapXlim, ErrorMapYlim, ErrorMapResolution)    
    ErrorMapL = np.full((ErrorMapResolution, ErrorMapResolution), np.finfo(np.float64).max, dtype=np.float64)
    for CamId in range(NumCams): ErrorMapL = np.minimum(ErrorMapL, ErrorMapsPolyL[CamId])
    
    fig = plt.figure(figsize=DefaultFigSize)
    ax  = fig.add_subplot(1,1,1)
    ax.set_title("Error map")
    a = VisualiserL.drawErrorMap(ax, ErrorMapL, ErrorMapXlim, ErrorMapYlim)
    VisualiserL.drawSystem(ax)
    ax.set(xlim=Xlim, ylim=Ylim)
    ax.grid()
    fig.tight_layout()    
    fig.savefig("LinearErr.pdf")
    if DisplayFigs: plt.show()

  if GenCir:  
    ErrorMapsPolyC = EstimatorC.calculateDeltaDistanceMapsPolygon(ErrorMapXlim, ErrorMapYlim, ErrorMapResolution) 
    ErrorMapC = np.full((ErrorMapResolution, ErrorMapResolution), np.finfo(np.float64).max, dtype=np.float64)
    for CamId in range(NumCams): ErrorMapC = np.minimum(ErrorMapC, ErrorMapsPolyC[CamId])

    fig = plt.figure(figsize=DefaultFigSize)
    ax  = fig.add_subplot(1,1,1)
    ax.set_title("Error map")
    VisualiserC.drawErrorMap(ax, ErrorMapC, ErrorMapXlim, ErrorMapYlim)
    VisualiserC.drawSystem(ax)
    ax.set(xlim=Xlim, ylim=Ylim)
    ax.grid()
    fig.tight_layout()
    fig.savefig("CircularErr.pdf")
    if DisplayFigs: plt.show()

if DrawCamMap:
  print("DrawErrorMap...")

  ErrorMapXlim = (-Distance, Distance)
  ErrorMapYlim = (- 1.2 * Distance, 1.2*Distance)

  print("calculateErrorMap...")  
  
  if GenLin:
    (ErrorMapsPolyL, CamMapsPolyL) = EstimatorL.calculateDeltaDistanceMapsSimplifiedForLinear(ErrorMapXlim, ErrorMapYlim, ErrorMapResolution)
    CamMapL = CamMapsPolyL[MainCamId]
    fig = plt.figure(figsize=DefaultFigSize)
    ax  = fig.add_subplot(1,1,1)
    ax.set_title("Best pair for main camera")
    VisualiserL.drawCamMap(ax, CamMapL, ErrorMapXlim, ErrorMapYlim)
    VisualiserL.drawSystem(ax)
    ax.set(xlim=Xlim, ylim=Ylim)
    ax.grid()
    plt.tight_layout()
    plt.savefig("LinearCam.pdf")
    if DisplayFigs: plt.show()

    CamMapLG = EstimatorL.recalculateDeltaDistanceMapsSimplifiedForLinearAndAllCameras(ErrorMapXlim, ErrorMapYlim, ErrorMapResolution,CamMapsPolyL)

    fig = plt.figure(figsize=DefaultFigSize)
    ax  = fig.add_subplot(1,1,1)
    ax.set_title("Number of baseline distances between camera pair")
    #VisualiserL.drawErrorMap(ax, CamMapLG, ErrorMapXlim, ErrorMapYlim)
    VisualiserL.drawCamMap(ax, CamMapLG, ErrorMapXlim, ErrorMapYlim, True)
    VisualiserL.drawCamDistances(ax,ErrorMapResolution)
    VisualiserL.drawSystem(ax)
    ax.set(xlim=Xlim, ylim=Ylim)
    ax.grid()
    plt.tight_layout()
    plt.savefig("LinearCamAll.pdf")
    if DisplayFigs: plt.show()

  if GenCir:

    (ErrorMapsPolyC, CamMapsPolyC) = EstimatorC.calculateDeltaDistanceMapsSimplifiedForCircular(ErrorMapXlim, ErrorMapYlim, ErrorMapResolution)
    CamMapC = CamMapsPolyC[MainCamId]

    fig = plt.figure(figsize=DefaultFigSize)
    ax  = fig.add_subplot(1,1,1)
    ax.set_title("Best pair for main camera")
    VisualiserC.drawCamMap(ax, CamMapC, ErrorMapXlim, ErrorMapYlim)
    VisualiserC.drawSystem(ax)
    ax.set(xlim=Xlim, ylim=Ylim)
    ax.grid()
    plt.tight_layout()
    plt.savefig("CircularCam.pdf")
    if DisplayFigs: plt.show()



if DrawSimplifiedErrLinear:
  print("DrawSimplifiedErrLinear...")

  ErrorMapXlim = (-Distance, Distance)
  ErrorMapYlim = (- 1.2 * Distance, 1.2*Distance)

  (ErrorMapsSimpleL, CamMapsSimpleL) = EstimatorL.calculateDeltaDistanceMapsSimplifiedForLinear(ErrorMapXlim, ErrorMapYlim, ErrorMapResolution)
  ErrorMapSimpleL = np.full((ErrorMapResolution, ErrorMapResolution), np.finfo(np.float64).max, dtype=np.float64)
  for ErrorMap in ErrorMapsSimpleL: 
    ErrorMapSimpleL = np.minimum(ErrorMapSimpleL, ErrorMap)

  fig = plt.figure(figsize=DefaultFigSize)
  ax  = fig.add_subplot(1,1,1)
  ax.set_title("Simplified error map")

  VisualiserL.drawErrorMap(ax, ErrorMapSimpleL, ErrorMapXlim, ErrorMapYlim)
  VisualiserL.drawSystem(ax)
  ax.set(xlim=Xlim, ylim=Ylim)
  ax.grid()

  plt.tight_layout()
  plt.savefig("LinearErrSimplified.pdf")
  if DisplayFigs: plt.show()

if DrawSimplifiedErrCircular:
  print("DrawSimplifiedErrCircular...")
  
  ErrorMapXlim = (-Distance, Distance)
  ErrorMapYlim = (- 1.2 * Distance, 1.2*Distance)

  (ErrorMapsSimpleC, CamMapsPolyC) = EstimatorC.calculateDeltaDistanceMapsSimplifiedForCircular(ErrorMapXlim, ErrorMapYlim, ErrorMapResolution)
  ErrorMapSimpleC = np.full((ErrorMapResolution, ErrorMapResolution), np.finfo(np.float64).max, dtype=np.float64)

  for ErrorMap in ErrorMapsSimpleC: 
    ErrorMapSimpleC = np.minimum(ErrorMapSimpleC, ErrorMap)

  fig = plt.figure(figsize=DefaultFigSize)
  ax  = fig.add_subplot(1,1,1)
  ax.set_title("Simplified error map")
  VisualiserL.drawErrorMap(ax, ErrorMapSimpleC, ErrorMapXlim, ErrorMapYlim)
  VisualiserC.drawSystem(ax)
  ax.set(xlim=Xlim, ylim=Ylim)
  ax.grid()
  plt.tight_layout()
  plt.savefig("CircularErrSimplified.pdf")
  if DisplayFigs: plt.show()
