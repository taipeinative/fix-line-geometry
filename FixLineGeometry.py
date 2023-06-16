import fiona
import geopandas as gpd
import pandas as pd
from shapely import LineString, MultiPoint, ops
from tqdm import tqdm

class FixLineGeometry:

    def __init__(self: 'FixLineGeometry', path: str, crs: int = 3826) -> None:
        '''
        The initial settings of FixLineGeometry.

        Parameters
        --------
        path: str
            The file path of the shapefile.

        crs: int
            Targeted crs for the following modifications. 
        '''
        geom = fiona.open(path)

        if (geom.schema['geometry'] == 'LineString'):   # Ref: https://gis.stackexchange.com/questions/457106/check-geometry-type-of-shapefile-with-geopandas

            self.NETWORK = gpd.read_file(path, encoding = 'utf-8').to_crs(crs)

        else:

            raise ValueError(f'Invalid geometry type: expect LineString but get {geom.schema["geometry"]}')

    def debug(
            
            self: 'FixLineGeometry',
            auto_fix: bool = True,
            threshold: float = 1.0,
            savefile_path: str = r'.',
            savefile_name: str = 'fixed'

        ) -> gpd.GeoDataFrame:
        '''
        List all not properly connected polylines.

        Parameters
        --------
        auto_fix: bool
            Whether to fix unconnected vetices or not. Default to be `False`.
        
        threshold: float
            The threshold to check whether the lines should be connected or not.
            For example, two non-connected vortex within 1 meters with threshold set to 1 meters would be regarded as unconnected.

        savefile_path: str
            The folder corrected shapefile should be placed at.

        savefile_name: str
            The name of the corrected shapefile. 
        '''

        def to_MultiPoint(feature: LineString) -> MultiPoint:
            '''
            Convert LineString into MultiPoint object; only apply to one row so use it with `lambda`.

            Parameters
            --------
            feature: shapely.LineString
                The linestring object needs converting into multipoint.

            Returns
            --------
            result: shapely.MultiPoint
                The multipoint version of original feature.
            '''

            return MultiPoint(feature.coords)

        multi_point_NETWORK = self.NETWORK.copy()
        # Convert the linestring into multipoint. / 將線轉換為複數點
        multi_point_NETWORK['geometry'] = self.NETWORK.apply(lambda L: to_MultiPoint(L['geometry']), axis = 1)

        problem_case = gpd.GeoDataFrame()

        # For auto fix mode / 自動修復模式用
        if (auto_fix):
            
            # Explode the multipoint so we can find out the point where the distance of any two line is smaller then threshold. /
            # 將複數點轉換成單一點，這樣才能個別檢查在哪一點任意兩條線的距離小於閥值。
            single_point_series = multi_point_NETWORK['geometry'].explode(index_parts = True).copy()
            # Perform buffer on every point with the radius using threshold. / 對每一點執行以閥值為半徑的環域。
            single_point_buffer = single_point_series.buffer(threshold)

            for (i, item) in tqdm(single_point_buffer.items(), 'Traversing every vertices', len(single_point_series)):

                # Find the linestrings which intersect the point buffer, or the linestring within the distance of threshold near the point. /
                # 找出和點環域有交集的線，即與點的直線距離在閥值內的線。
                close_NETWORK = self.NETWORK[self.NETWORK['geometry'].intersects(item)].copy()
                # Find the line itself. Since the single points were converted from multipoints, which were extracted from the linestrings, /
                # there will always be the original line in the result. /
                # 找出線本身。由於每一單點是由複數點轉換而來，而複數點是線上的每一頂點，因此必定會搜尋到原始的線。
                close_NETWORK['self'] = close_NETWORK['geometry'].geom_equals(self.NETWORK.at[i[0], 'geometry'])
                # Find the linestings touching but not crossing the original line. / 找出碰觸但未超出原始線的線。
                close_NETWORK['touch'] = close_NETWORK['geometry'].touches(self.NETWORK.at[i[0], 'geometry'])

                def check_within1(idx: tuple, feature: LineString) -> bool:
                    '''
                    Check the current point is within the linestring; only apply to one row so use it with `lambda`.

                    Parameters
                    --------
                    idx: tuple
                        The current index of the base point in `single_point_series` (the point we want to check within).

                    feature: shapely.LineString
                        The linestring to compare with the base point.

                    Returns
                    --------
                    result: bool
                        Whether the point is within the linestring.
                    '''

                    base_point = single_point_series[idx[0]][idx[1]]
                    ref_multipoint = to_MultiPoint(feature)

                    return base_point.within(ref_multipoint)

                # Find out whether the intersection point of any two linestrings are correctly placed. /
                # Note that this only catch the point between the start and finish of the line; /
                # Remaing case have been catched by previous `touch`. /
                # 確認任意兩線交會處的頂點有正確的疊合。這只適用在線的端點外的情況，端點的連接情況已在`touch`欄位確認。
                close_NETWORK['within'] = close_NETWORK.apply(lambda L: check_within1(i, L['geometry']), axis = 1)

                # Exclude the record if any of these columns are `True`: self, touch & within. /
                # This means the record is connect correctly with the line which the point buffer belongs to. /
                # 排除以下任意為真的欄位：self、touch 和 within，這代表這筆紀錄的線與點環域所屬的線有正確連接上。
                close_NETWORK = close_NETWORK.query('(self == False) & (touch == False) & (within == False)').copy()
                # Drop the column we don't need anymore. / 刪除接下來用不到的欄位。
                close_NETWORK.drop(['self', 'touch', 'within'], axis = 1, inplace = True)

                # Record problem pairs / 記錄未正確連接的線組合。
                if (len(close_NETWORK) > 0):

                    geom_idx = close_NETWORK.columns.get_loc('geometry')
                    close_NETWORK.insert(loc = geom_idx,     column = 'nbr_idx0', value = i[0])
                    close_NETWORK.insert(loc = geom_idx + 1, column = 'nbr_idx1', value = i[1])
                    close_NETWORK.insert(loc = geom_idx + 2, column = 'nbr_name', value = self.NETWORK.at[i[0], 'name'])

                    if ('section' in self.NETWORK.columns):

                        close_NETWORK.insert(loc = geom_idx + 3, column = 'nbr_sctn', value = self.NETWORK.at[i[0], 'section'])
                        close_NETWORK.insert(loc = geom_idx + 4, column = 'nbr_poin', value = single_point_series[i[0]][i[1]])

                    else:

                        close_NETWORK.insert(loc = geom_idx + 3, column = 'nbr_poin', value = single_point_series[i[0]][i[1]])

                    problem_case = pd.concat([problem_case, close_NETWORK]).copy()

        # For non auto fix mode; basically same as above but use a more simple calculaing method. /
        # CANNOT EXCLUDE GENERAL INTERSECTIONS! USE AUTO FIX MODE FOR THE PURPOSE.
        # 非自動修復模式用；基本上和自動修復模式相同，但使用較簡單的算法。無法排除兩條交會點在線段中間的點，即使正確也將視為錯誤。
        else:

            multi_point_buffer = multi_point_NETWORK.copy()
            multi_point_buffer['geometry'] = multi_point_NETWORK['geometry'].buffer(threshold)

            for (i, row) in tqdm(multi_point_buffer.iterrows(), 'Traversing every line', len(multi_point_buffer)):

                close_NETWORK = self.NETWORK[self.NETWORK['geometry'].intersects(row['geometry'])].copy()
                close_NETWORK['self'] = close_NETWORK['geometry'].geom_equals(self.NETWORK.at[i, 'geometry'])
                close_NETWORK['touch'] = close_NETWORK['geometry'].touches(self.NETWORK.at[i, 'geometry'])

                close_NETWORK = close_NETWORK.query('(touch == False) & (self == False)').copy()
                close_NETWORK.drop(['self', 'touch'], axis = 1, inplace = True)

                if (len(close_NETWORK) > 0):

                    geom_idx = close_NETWORK.columns.get_loc('geometry')
                    close_NETWORK.insert(loc = geom_idx,     column = 'nbr_idx0', value = i)
                    close_NETWORK.insert(loc = geom_idx + 1, column = 'nbr_name', value = self.NETWORK.at[i, 'name'])

                    if ('section' in self.NETWORK.columns):

                        close_NETWORK.insert(loc = geom_idx + 2, column = 'nbr_sctn', value = self.NETWORK.at[i, 'section'])

                    problem_case = pd.concat([problem_case, close_NETWORK]).copy()

        # Append new index to all problem cases. / 為所有錯誤的線組合添加新的索引。
        problem_case.reset_index(inplace = True)

        problem_NETWORK = gpd.GeoDataFrame()
        duplicate_filter = set()

        # Remove the duplicates. Since the function iterrated all point buffer, any incorrectly connected network would append new pairs to the `problem_case`. /
        # If there are two incorrectly connected lines, there will be two records where one of the item is a reversed verion of another one. /
        # 移除重複的線組合。由於上方的方式迭代了每一個點環域，任何未正確連接的線皆會產生一次記錄，導致記錄的數量是實際問題的兩倍；/
        # 舉例來說，甲線的A點和乙線的B點距離在閥值內，迭代完成後將會出現 (甲, 乙) 和 (乙, 甲) 的組合，其中一組會是另一組的相反。 
        for (i, row) in tqdm(problem_case.iterrows(), 'Removing duplicate networks', len(problem_case)):

            item = (row['index'], row['nbr_idx0'])

            # If the record doesn't exist in the `duplicate_filter` yet, add it and its reversed version so next time its reversed version will be omitted. /
            # 若該筆記錄尚未存在於 `duplicate_filter`中，將其與其相反版本加入，如此一來其相反版本將不會被加入。
            if ((item not in duplicate_filter) & (item[::-1] not in duplicate_filter)):

                duplicate_filter.add(item)
                duplicate_filter.add(item[::-1])

                # ... so that `problem_NETWORK` will not have duplictes. /
                # ……因此 `problem_NETWORK` 中不會有重複的紀錄。
                problem_NETWORK = pd.concat([problem_NETWORK, problem_case.loc[[i]]]).copy()

        if (not auto_fix):

            return problem_NETWORK

        else:

            fixed_point_series = single_point_series.copy()

            # Find the nearest point of the two lines and change the coordinates. / 找出兩條線相交的最近點並將其校正。
            for (i, row) in tqdm(problem_NETWORK.iterrows(), 'Searching nearest vertices', len(problem_NETWORK)):

                host_point = row['nbr_poin']
                guest_line = multi_point_NETWORK.at[row['index'], 'geometry']

                (query, nearest) = ops.nearest_points(host_point, guest_line)

                guest_point = single_point_series[single_point_series == nearest]
                guest_index = guest_point.index

                fixed_point_series[guest_index] = host_point

            # Flatten the multi-index series. Since it was originally exploded from multipoint, /
            # it is related to the linestring in original NETWORK dataframe./
            # 移除 series 的複數索引。由於其複數索引來自多點，因此新的索引和NETWORK中的線的索引相連。
            fixed_point_series = fixed_point_series.droplevel(1, axis = 0)

            fixed_point_list = list()

            # Implode the series so the series is made of multipoint. / 重新將單點組裝回多點。
            for (i, item) in fixed_point_series.items():

                if (len(fixed_point_list) == i):

                    fixed_point_list.append([item])

                else:

                    fixed_point_list[i].append(item)

            fixed_line_list = list()

            # Convert multipoint into linestring / 將多點轉換成線。
            for multi_point in fixed_point_list:

                fixed_line_list.append(LineString(multi_point))

            # Replace old NETWORK's geometry. / 替換掉舊 NETWORK 中的幾何。
            fixed_NETWORK = self.NETWORK.copy()
            fixed_NETWORK['geometry'] = fixed_line_list

            # Export to shapefile. / 輸出shapefile檔案。
            fixed_NETWORK.to_file(f'{savefile_path}\\{savefile_name}.shp', encoding = 'utf-8')

            return problem_NETWORK

    def plot(self: 'FixLineGeometry'):
        '''
        Plot the map of the shapefile.
        '''
        basemap = self.NETWORK.plot(color = 'lightgrey', linewidth = 0.75)
        basemap.set_axis_off()
        basemap.set_title('Network Map')

    # Properties

    @property
    def NETWORK(self: 'FixLineGeometry') -> gpd.GeoDataFrame:
        '''
        The network GeoDataFrame of the file.
        '''
        return self._NETWORK
    
    @NETWORK.setter
    def NETWORK(self: 'FixLineGeometry', value):
        '''
        Set the network GeoDataFrame.
        '''
        self._NETWORK = value