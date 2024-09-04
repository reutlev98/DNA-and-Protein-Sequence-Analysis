# Reut Lev 207385741
import pandas as pd


class myData:
    def __init__(self, path) -> None:
        self.data = pd.read_csv(path, low_memory=False)

    # EX1
    def num_player_height(self, x, y):
        """
        Given height x & y (x<y) the function returns the number of players whose height is between x & y, inclusive.
        :param x: height (x<y)
        :param y: height (x<y)
        :return: number of players whose height is between x & y, inclusive
        """
        filter_data = self.data[self.data["height_cm"].between(x, y)]  # filter data where height col is between x & y
        number_of_player = filter_data.shape[0]  # counting the number of rows
        return number_of_player

    # EX2
    def df_birthyear(self, year):
        """
        the function returns a dataframe with the name_short, name_club columns of the players born in the given year.
        :param year: The year the actor was born
        :return:dataframe with the name_short, name_club columns of the players born in the given year
        """
        self.data["dob"] = pd.to_datetime(self.data["dob"], format="%Y/%m/%d")  # parse dob column into datetime objects
        self.data["year_of_birth"] = self.data["dob"].dt.year  # uses dt to extract the year from each datetime object
        data_filterd = self.data[self.data["year_of_birth"] == year][["short_name", "club_name"]]  # filter given year
        return data_filterd

    # EX3
    def list_sorted(self, col1, col2, k):
        """
        :param col1: first column
        :param col2: second column
        :param k: a number
        :return: list of the k players with the highest values in the first given column.The list will include
        the name_short of the players. The list will be sorted by the values in the second column given in ascending
        order from the lowest value to the highest value.
        """

        data_sorted = self.data.sort_values([col1, col2], ascending=[False, True])  # sorted desc order by col1 values
        temp_data = data_sorted.head(k)  # Filtering the information according to the highest K's
        data_sorted = temp_data.sort_values([col2], ascending=[True])  # sorted asce order by col1 values
        player_list = data_sorted["short_name"].head(k).tolist()  # convert the resulting series to a list
        return player_list

    # EX4
    def tuples_players_by_year(self, x, y):
        """
        :param x: number of year (x<y)
        :param y: number of year (x<y)
        :return: a list of tuples of the year and the number of players born in that year,starting from year X to Y.
        """
        self.data["dob"] = pd.to_datetime(self.data["dob"], format="%Y/%m/%d")  # parse dob column into datetime objects
        self.data["dob"] = self.data["dob"].dt.year  # uses dt to extract the year from each datetime object
        year_counts = self.data["dob"].value_counts().sort_index()  # sorts the resulting series by the index (year)
        year_list = [(year, count) for year, count in year_counts.items() if x <= year <= y]  # create a list of tuples
        return year_list

    # EX5
    def mean_std(self, col, first_name):
        """
        :param col: a column
        :param first_name: first name
        :return: the mean and standard deviation,of the values in the given column for the players with the first name.
        """
        self.data["first_name"] = self.data["long_name"].str.split(" ").str[0]  # split the name column, selects str[0]
        filtered_data = self.data[self.data["first_name"] == first_name][[col]]  # select only the rows value=first name
        mean = filtered_data[col].mean()  # calculates the mean
        std = filtered_data[col].std()  # calculates standard deviation
        return mean, std

    # EX6
    def max_players(self, col):
        """
        :param col: a column
        :return: the value in this column with the largest number of players (the value with the most repetitions).
        """
        value_count = self.data[col].value_counts()  # count the number of occurrences of each value in the col column
        max_value = value_count.idxmax()  # find the value with the highest count
        return max_value
