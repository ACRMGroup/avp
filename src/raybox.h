BOOL DoesLineCrossBox(VEC3F start, VEC3F stop, REAL length, 
                      VEC3F corner1, VEC3F corner2);
BOOL HitBoundingBox(REAL minB[3], REAL maxB[3],
                    REAL origin[3], REAL dir[3],
                    REAL coord[3]);
